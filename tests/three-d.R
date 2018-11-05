rot_x <- function(deg) {
  r <- deg / 180 * pi
  matrix(byrow=TRUE, nrow=3, c(
         1,       0,       0,
         0,  cos(r), -sin(r),
         0,  sin(r),  cos(r))
  )
}
rot_y <- function(deg) {
  r <- deg / 180 * pi
  matrix(byrow=TRUE, nrow=3, c(
    cos(r),       0,  sin(r),
         0,       1,       0,
   -sin(r),       0,  cos(r))
  )
}
rot_z <- function(deg) {
  r <- deg / 180 * pi
  matrix(byrow=TRUE, nrow=3, c(
    cos(r), -sin(r),       0,
    sin(r),  cos(r),       0,
         0,       0,       1)
  )
}

# eltif <- raster::raster("~/Downloads/dem_01.tif")
# elmat1 <- matrix(
#   raster::extract(eltif,raster::extent(eltif),buffer=10000),
#   nrow=ncol(eltif),ncol=nrow(eltif)
# )
# mx2 <- volcano
# mx2 <- elmat1
# sun <- 135
# els <- seq(-90, 90, length=25)
# sh2 <- rayshader::ray_shade(mx2, els, sun, lambert=FALSE)
# sh2 <- sh2[, rev(seq_len(ncol(sh2)))]

# Convert matrix to long format, and rotate

mxl <- cbind(x=c(row(mx2)), y=c(col(mx2)), z=c(mx2))
mxlr <- mxl %*% t(rot_z(215)) %*% t(rot_x(-30))

# Account for perspective; this is not correct as we're doing it purely on the
# basis of how far from the observer we are along the y axis, instead of actual
# distance.

obs <- c(
  diff(range(mxlr[,1])) / 2 + min(mxlr[,1]),
  diff(range(mxlr[,2])) / 2 + min(mxlr[,2]),
  diff(range(mxlr[,3])) + max(mxlr[,3])
)
dist.to.obs <- sqrt(colSums((t(mxlr) - obs)^2))

fovh <- 60 / 180 * pi
fovv <- 60 / 180 * pi

mxlrp <- mxlr
# mxlrp[,1] <- mxlrp[,1] / ((mxlrp[,2] - obs[,2]) * tan(fovh) / 2)
# mxlrp[,3] <- mxlrp[,3] / ((mxlrp[,2] - obs[,2]) * tan(fovv) / 2)

# ggplot(as.data.frame(mxlrp)) + geom_point(aes(V1, V2, color=c(sh2)))

# Project onto a canvas, we will do this by interpreting our grid of points as a
# set of tiles where the points represent the vertices of each tile.  Start by
# transforming our point data into tile data

res.x <- 600
res.y <- 500

# Add meta data, and rescale x,y values into 1:res.x and 1:res.y

mxres <- cbind(mxlrp, shadow=c(sh2), dist=dist.to.obs)
mxres[,2] <- (mxres[,2] - min(mxres[,2])) / diff(range(mxres[,2])) *
  (res.y - 1) + 1
mxres[,1] <- (mxres[,1] - min(mxres[,1])) / diff(range(mxres[,1])) *
  (res.x - 1) + 1

# Generate coordinates for triangular mesh from our points.

triangular_mesh <- function(dat, nr, nc) {
  stopifnot(identical(dim(dat)[2], 5L), nr * nc == nrow(dat))

  # create a tile mesh by combining four copies of our original matrix,
  # dropping first/last row/col in turn.  The arithmetic is to map location
  # of first/last row/col to our long input matrix `dat`.

  lr.lc <- -c(seq_len(nc) * nr, seq_len(nr) + nr * (nc - 1L))
  lr.fc <- -c(seq_len(nc) * nr, seq_len(nr))
  fr.lc <- -c((seq_len(nc) - 1L) * nr + 1L, seq_len(nr) + nr * (nc - 1L))
  fr.fc <- -c((seq_len(nc) - 1L) * nr + 1L, seq_len(nr))

  tile.mesh <- c(
    dat[lr.lc,],   # drop last row, last col
    dat[lr.fc,],   # drop last row, first col
    dat[fr.lc,],   # drop first row, last col
    dat[fr.fc,]    # drop first row, first col
  )
  # store mesh data in a 3D array, where the indices of the 3rd dimension
  # represent which vertex of the mesh element is described.

  dim(tile.mesh) <- c(length(tile.mesh) / (4L * ncol(dat)), ncol(dat), 4L)

  # now split tiles into triangles

  mesh <- array(0, c(nrow(tile.mesh)*2,ncol(tile.mesh),3))
  mesh[seq_len(nrow(tile.mesh)),,] <- tile.mesh[,,1:3]
  mesh[-seq_len(nrow(tile.mesh)),,] <- tile.mesh[,,2:4]
  mesh
}
mesh <- triangular_mesh(mxres, nr=nrow(mx2), nc=ncol(mx2))

# For each triangle in the mesh, determine which points in the integer grid
# could potentially be in the mesh

candidate_points <- function(mesh) {
  stopifnot(identical(dim(mesh)[2:3], c(5L,3L)))

  # assume all points between max and min on y/x could be in mesh
  points.int <- cbind(
    id=seq_len(nrow(mesh)),
    xmin=pmin(mesh[,1L,1L], mesh[,1L,2L], mesh[,1L,3L]),
    xmax=pmax(mesh[,1L,1L], mesh[,1L,2L], mesh[,1L,3L]),
    ymin=pmin(mesh[,2L,1L], mesh[,2L,2L], mesh[,2L,3L]),
    ymax=pmax(mesh[,2L,1L], mesh[,2L,2L], mesh[,2L,3L])
  )
  # don't allow the maxs to be exactly at the boundary

  xmax.bad <- points.int[,'xmax'] == floor(points.int[,'xmax'])
  points.int[xmax.bad, 'xmax'] <- points.int[xmax.bad, 'xmax'] -
    (points.int[xmax.bad, 'xmax'] - points.int[xmax.bad, 'xmin']) * .01
  ymax.bad <- points.int[,'ymax'] == floor(points.int[,'ymax'])
  points.int[ymax.bad, 'ymax'] <- points.int[ymax.bad, 'ymax'] -
    (points.int[ymax.bad, 'ymax'] - points.int[ymax.bad, 'ymin']) * .01

  points.int[, c('xmin', 'ymin')] <- ceiling(points.int[, c('xmin', 'ymin')])
  points.int[, c('xmax', 'ymax')] <- floor(points.int[, c('xmax', 'ymax')])

  # for each triangle, generate the set of integer coordinates that it could
  # contain.

  p.int.len.x <- pmax(points.int[, 'xmax'] - points.int[, 'xmin'] + 1, 0)
  p.int.len.y <- pmax(points.int[, 'ymax'] - points.int[, 'ymin'] + 1, 0)
  p.int.len <- p.int.len.x * p.int.len.y

  # For each id, need the equivalent of `expand.grid` of the xmin:xmax and
  # ymin:ymax sequences.  The following logic replicates that by doing it for
  # 1:(xmax-xmin) and 1:(ymax-ymin); we add back xmin and xmax later.  Contorted
  # logic is to minimize calls to R functions as this is bottleneck.

  seq.id <- rep(points.int[, 'id'], p.int.len)
  seq.x <- rep(sequence(p.int.len.x), rep(p.int.len.y, p.int.len.x)) - 1L
  seq.x.delta <- which(c(FALSE, diff(seq.x) != 0L | diff(seq.id) != 0L))

  seq.y <- integer(length(seq.x))
  seq.y[seq.x.delta] <- seq.x.delta - 1L
  seq.y <- seq_along(seq.y) - cummax(seq.y)

  # add back xmin and xmax

  cbind(
    id=seq.id,
    x.int=seq.x + points.int[seq.id, 'xmin'],
    y.int=seq.y + points.int[seq.id, 'ymin']
  )
}
points.raw <- candidate_points(mesh)
#
# Remove all points that are not actually within their mesh triangle

trim_points <- function(points, mesh) {
  stopifnot(ncol(points) == 3, identical(dim(mesh)[-1], c(5L,3L)))

  # For each point, compute the vectors from that point to each of the vertices
  # of the mesh triangle we are checking it is in

  vecs <- mesh[points[,'id'],1:2,] - array(points[,2:3], c(nrow(points),2,3))

  # compute angles between the triangle vertices: a . b = |a||b| cos(p)

  sqrRsq <- function(x) sqrt(rowSums(x^2))

  mesh.ang.1 <-
    acos(rowSums(vecs[,,1]*vecs[,,2]) / (sqrRsq(vecs[,,1])*sqrRsq(vecs[,,2])))
  mesh.ang.2 <-
    acos(rowSums(vecs[,,2]*vecs[,,3]) / (sqrRsq(vecs[,,2])*sqrRsq(vecs[,,3])))
  mesh.ang.3 <-
    acos(rowSums(vecs[,,3]*vecs[,,1]) / (sqrRsq(vecs[,,3])*sqrRsq(vecs[,,1])))

  # points that are within their triangles will have a sum of angles equal to
  # 2pi Need to be a little more generous for precision pruposes, could lead to
  # incorrect determinations in very corner cases

  points[(mesh.ang.1 + mesh.ang.2 + mesh.ang.3 >= 2 * pi - 1e-6),]
}
points.t <- trim_points(points.raw, mesh)

# Compute meta data for each point
#
# We need the shadow value, as well as the distance from observer.  We will do
# the distance weighted average of the mesh triangle vertex values which contain
# this information.

points_meta <- function(p, m) {
  m.dist <- sqrt(
    rowSums(
      aperm(
        (m[p[,'id'],1:2,] - array(p[,2:3],c(nrow(p),2,3))) ^ 2,
        c(1, 3, 2)
      ),
      dims=2
  ) )
  p.shadow <- rowSums(m[p[,'id'], 4,] * (1/m.dist)) / rowSums(1/m.dist)
  p.dist <- rowSums(m[p[,'id'], 5,] * (1/m.dist)) / rowSums(1/m.dist)
  cbind(p, shadow=p.shadow, dist=p.dist)
}
points.dat <- points_meta(points.t, mesh)

# Drop points hidden by others
#
# We use distance from observer, and keep only closest points for any given
# x,y coordinate.  Basically, we order all points by distance, and then drop any
# duplicated x,y coordinates so that we just keep the first of each.

drop_hidden_points <- function(p) {
  in.id <- seq_len(nrow(p))
  dist.ord <- order(p[,'dist'])

  id.unique <- -in.id[dist.ord][duplicated(p[dist.ord, 2:3])]
  p[id.unique,,drop=FALSE]
}
points <- drop_hidden_points(points.dat)

# mesh2 <- aperm(mesh, c(1, 3, 2))
# mesh3 <- cbind(
#   rep(seq_len(nrow(mesh2)), each=ncol(mesh2)),
#   c(t(mesh2[,,1])), c(t(mesh2[,,2]))
# )
# Need to reorder points in the mesh

# ggplot(as.data.frame(mesh3), aes(V2, V3, group=V1)) +
#  geom_polygon(color='grey', size=.2)

# compute vectors from integer points to the triangle vertices

# Drop duplicates x-y values; we keep those closest to the observer; use
# untransformed

stop()
ggplot(data=as.data.frame(points)) +
  geom_point(aes(x=x.int, y=y.int, group=NULL, color=shadow)) +
  scale_color_gradient(low='#333333', high='#ffffff', guide=FALSE)

#ggplot(as.data.frame(mesh3), aes(x=V2, y=V3, group=V1)) +
  # geom_polygon(size=.2, alpha=.4) +
  # scale_color_gradient(low='#333333', high='#ffffff', guide=FALSE)
  # coord_cartesian(xlim=c(26,30), ylim=c(15,20))
#   geom_polygon(
#     data=as.data.frame(t(mesh[9718,,][1:2,])), aes(V1, V2, group=21435145123),
#     color='red', fill='blue'
#   )


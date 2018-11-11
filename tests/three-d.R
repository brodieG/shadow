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
triangular_mesh <- function(dat, nr, nc) {
  stopifnot(identical(dim(dat)[2], 4L), nr * nc == nrow(dat))

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
  dimnames(mesh) <- list(NULL, c('x', 'y', 'z', 'texture'))
  mesh
}
# For each triangle in the mesh, determine which points in the integer grid
# could potentially be in the mesh

candidate_points <- function(mesh) {
  stopifnot(identical(dim(mesh)[2:3], c(4L,3L)))

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
  seq.y <- seq_along(seq.y) - cummax(seq.y) - 1L

  # add back xmin and xmax

  cbind(
    id=seq.id,
    x=seq.x + points.int[seq.id, 'xmin'],
    y=seq.y + points.int[seq.id, 'ymin']
  )
}
#
# Remove all points that are not actually within their mesh triangle

trim_points <- function(points, mesh) {
  stopifnot(ncol(points) == 3, identical(dim(mesh)[-1], c(4L,3L)))

  # For each point, compute the vectors from that point to each of the vertices
  # of the mesh triangle we are checking it is in

  vecs <- mesh[points[,'id'],1:2,] - array(points[,2:3], c(nrow(points),2,3))

  # compute angles between the triangle vertices: a . b = |a||b| cos(p)

  sqrRsq <- function(x) sqrt(rowSums(x*x))
  v1 <- vecs[,,1]
  v2 <- vecs[,,2]
  v3 <- vecs[,,3]
  v1s <- sqrRsq(v1)
  v2s <- sqrRsq(v2)
  v3s <- sqrRsq(v3)

  mesh.ang.1 <- acos(rowSums(v1*v2) / (v1s*v2s))
  mesh.ang.2 <- acos(rowSums(v2*v3) / (v2s*v3s))
  mesh.ang.3 <- acos(rowSums(v3*v1) / (v3s*v1s))

  # points that are within their triangles will have a sum of angles equal to
  # 2pi Need to be a little more generous for precision pruposes, could lead to
  # incorrect determinations in very corner cases

  points[(mesh.ang.1 + mesh.ang.2 + mesh.ang.3 >= 2 * pi - 1e-6),]
}
#' Compute Barycentric Coordinates
#'
#' Taken from [wikipedia](https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Conversion_between_barycentric_and_Cartesian_coordinates).
#'
#' @param p numeric matrix with two columns 'x' and 'y': the points to compute
#'   the barycentric coordinates for.
#' @param v numeric array with dim(n,2,3) where the columns are 'x' and 'y', and
#'   the 3rd dimension represents each vertex of the triangle.

barycentric_coord <- function(p, v) {
  x1 <- v[,'x',1]
  x2 <- v[,'x',2]
  x3 <- v[,'x',3]
  y1 <- v[,'y',1]
  y2 <- v[,'y',2]
  y3 <- v[,'y',3]
  x <- p[,'x']
  y <- p[,'y']

  y2.y3 <- y2 - y3
  x3.x2 <- x3 - x2
  x1.x3 <- x1 - x3
  x.x3 <- x - x3
  y.y3 <- y - y3
  det <- y2.y3 * x1.x3 + x3.x2 * (y1 - y3)

  l1 <- (y2.y3 * x.x3 + x3.x2 * y.y3) / det
  l2 <- ((y3 - y1) * x.x3 + x1.x3 * y.y3) / det
  l3 <- 1 - l1 - l2

  cbind(l1, l2, l3)
}
# Compute meta data for each point
#
# We need the shadow value, as well as the distance from observer.  We will do
# the distance weighted average of the mesh triangle vertex values which contain
# this information.

points_meta <- function(p, m) {
  # compute barycentric coordinates for each point relative to the vertices

  m.dat <- m[p[, 'id'], c('x', 'y', 'z', 'texture'),]
  bary <- barycentric_coord(p[, c('x', 'y')], m.dat)
  p.texture <- rowSums(bary * m.dat[,'texture',])
  p.z <- rowSums(bary * m.dat[,'z',])
  cbind(p, z=p.z, texture=p.texture)
}
empty_rows <- function(arr) (rowSums(arr) == -3 * ncol(arr))
empty_cols <- function(arr) (rowSums(colSums(arr)) == -3 * nrow(arr))

#' @param resolution integer(2L) output resolution where the first value is the
#'   x resolution, and the second the y resolution.
#' @param parallax numeric vector of degrees to add to the z angle of
#'   `rotation`.  A view will be rendered for each angle and returned in the
#'   result list.  This makes it easy to generate stereoscopic images.
#' @param dist numeric(1L) distance of the observer from the x-y plane as a
#'   multiple of the z distance spanned by the elevation surface after rotation.
#'   currently the observer is always centered on the x-y plane.
#' @param fov numeric(1L) field of view in degrees.

project_elev <- function(
  elevation, texture, rotation, parallax=0, dist=2,
  resolution=as.integer(dim(elevation) * 1.2)
) {
  mxl <- rbind(c(row(elevation)), c(col(elevation)), c(elevation))
  res <- vector('list', length(parallax))

  for(i in seq_along(parallax)) {
    mxlr <- rotation %*% rot_z(parallax[i]) %*% (mxl)
    ranges <- apply(mxlr, 1, range)
    obs <- c(
      diff(ranges[,1])/2 + ranges[1,1],
      diff(ranges[,2])/2 + ranges[1,2],
      dist * diff(ranges[,3]) + ranges[2,3]
    )
    mxlro <- mxlr - obs
    dimnames(mxlro)[[1]] <- c('x', 'y', 'z')

    # apply parallax rotation and compute observer position
    # compute projection

    mxlrpp <- mxlro
    mxlrpp['x',] <- -1/mxlro['z',] * mxlro['x',]
    mxlrpp['y',] <- -1/mxlro['z',] * mxlro['y',]

    # df <- as.data.frame(cbind(t(mxlrpp), color=c(sh2)))
    # print(ggplot(df[order(df$V3),]) + geom_point(aes(V1, V2, color=color)))
    # print(ggplot(df[order(df$V3),]) + geom_point(aes(V1, V2, color=V3)))
    # mxlrp[,3] <- mxlrp[,3] / ((mxlrp[,2] - obs[,2]) * tan(fovv) / 2)

    # ggplot(as.data.frame(t(mxlr))) + geom_point(aes(V1, V2, color=c(sh2)))
    # ggplot(as.data.frame(t(mxlrpp))) + geom_point(aes(V1, V2, color=c(sh2)))
    # ggplot(as.data.frame(cbind(t(mxl[1:2,]), z=c(texture)))) +
    #   geom_point(aes(x, y, color=z))

    # Add meta data, and rescale x,y values into 1:res.x and 1:res.y

    mxres <- cbind(t(mxlrpp), texture=c(texture))
    mxres[,'y'] <- (mxres[,'y'] - min(mxres[,'y'])) / diff(range(mxres[,'y'])) *
      (resolution[2] - 1) + 1
    mxres[,'x'] <- (mxres[,'x'] - min(mxres[,'x'])) / diff(range(mxres[,'x'])) *
      (resolution[1] - 1) + 1

    # print(ggplot(as.data.frame(mxres)) + geom_point(aes(V1, V2, color=texture)))

    # Generate coordinates for triangular mesh from our points.

    mesh <- triangular_mesh(mxres, nr=nrow(elevation), nc=ncol(elevation))

    # mesh2 <- aperm(mesh, c(1, 3, 2))
    # mesh3 <- cbind(
    #   rep(seq_len(nrow(mesh2)), each=ncol(mesh2)),
    #   c(t(mesh2[,,1])), c(t(mesh2[,,2])), c(t(mesh2[,,5]))
    # )
    # meshdf <- as.data.frame(mesh3)
    # meshdf$z <- ave(meshdf$V4, meshdf$V1, FUN=mean)
    # meshdf$id <- rep(seq_len(nrow(meshdf)/3), each=3)
    # meshdf$id <- factor(meshdf$id, levels=sort(unique(meshdf$id), dec=T))
    # ggplot(meshdf, aes(V2, V3, group=id)) +
    #  geom_polygon(color='grey', size=.05)

    points.raw <- candidate_points(mesh)
    points.t <- trim_points(points.raw, mesh)
    points.dat <- points_meta(points.t, mesh)
    # resolve duplicates by putting nearest points last
    points <- points.dat[order(points.dat[,'z']),]

    # ggplot(as.data.frame(points.dat)) +
    # geom_point(aes(x=x.int, y=y.int, color=texture))

    res.mx <- matrix(-1, nrow=resolution[1], ncol=resolution[2])
    res.mx[points[, c('x', 'y')]] <- points[, 'texture']

    # transformations so renders okay in png

    res.mx <- t(res.mx)
    res.mx <- res.mx[rev(seq_len(nrow(res.mx))), ]
    res[[i]] <- array(res.mx, c(dim(res.mx), 3L))
  }
  res
}


# eltif <- raster::raster("~/Downloads/dem_01.tif")
# elmat1 <- matrix(
#   raster::extract(eltif,raster::extent(eltif),buffer=10000),
#   nrow=ncol(eltif),ncol=nrow(eltif)
# )
mx2 <- volcano
# mx2 <- elmat1
sun <- -45
els <- seq(-90, 90, length=25)
sh2 <- rayshader::ray_shade(mx2, els, sun, lambert=FALSE)
sh2 <- sh2[, rev(seq_len(ncol(sh2)))]
sh2 <- sh2 * .9 + .1

# Convert matrix to long format, and rotate

angle <- 3 * c(1,-1)
#rot <- rot_x(-20) %*% rot_z(215)
rot <- rot_x(-20) %*% rot_z(65)
#proj <- project_elev(mx2, sh2, rot, parallax=angle)
stop()
treeprof(
proj <- project_elev(
  mx2, sh2, rot, parallax=angle, dist=.5, resolution=c(800,800)
)
)
system.time(
proj <- project_elev(mx2, sh2, rot, parallax=0, dist=.5, resolution=c(800,800))
)
left <- proj[[1]]
right <- proj[[2]]

empty.rows <- which(empty_rows(left) & empty_rows(right))
empty.cols <- which(empty_cols(left) & empty_cols(right))
left <- left[-empty.rows, -empty.cols, ]
png::writePNG(left, 'persp-left.png')
right <- right[-empty.rows, -empty.cols, ]

# determine what offset minimizes the mismatch of the max y values.

offset <- 0
margin <- 10

# right[right < 0] <- .8
right[,,1] <- 0
# left[left < 0] <- .8
left[,,2:3] <- 0

# right[,,1] <- 0

res <- array(0, dim(left) + c(margin*2, offset + margin*2, 0))
rows <- seq_len(nrow(left))
cols <- seq_len(ncol(left))

res[rows + margin, cols + margin, ] <- right
res[rows + margin, cols + offset + margin, ] <-
  pmin(res[rows + margin, cols + offset + margin, ] + left, 1)
png::writePNG(res, 'persp-color.png')

# side by side

proj <- project_elev(mx2, sh2, rot, parallax=angle, dist=.5, res=c(400,350))
left <- proj[[1]]
right <- proj[[2]]

empty.rows <- which(empty_rows(left) & empty_rows(right))
empty.cols <- which(empty_cols(left) & empty_cols(right))

# drop extra channels

left <- left[-empty.rows, -empty.cols, 1]
right <- right[-empty.rows, -empty.cols, 1]

no.overlap <- left[, rev(seq_len(ncol(left)))] == -3 | right == -3
no.over.rle <- rle(colSums(no.overlap) == nrow(no.overlap))
overlap.size <- if(isTRUE(no.over.rle[['values']][1]))
  no.over.rle[['lengths']][1] else 0L

# combine and write png

combined <- array(0, dim=dim(left) * c(1,2) - c(0, overlap.size))
combined[,seq_len(ncol(left))] <- left
combined[,tail(seq_len(ncol(combined)), ncol(right))] <-
  combined[,tail(seq_len(ncol(combined)), ncol(right))] + right

margin <- c(100, 20)
combined.fin <- array(0, dim(combined) + margin * 2)
combined.fin[
  seq_len(nrow(combined)) + margin[1],
  seq_len(ncol(combined)) + margin[2]
] <- combined
png::writePNG(combined.fin, 'persp.png')


stop('done')


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
library(ggplot2)
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


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


mx2 <- volcano
sun <- 135
els <- seq(-90, 90, length=25)
sh2 <- rayshader::ray_shade(mx2, els, sun, lambert=FALSE)
sh2 <- sh2[, rev(seq_len(ncol(sh2)))]

mxl <- cbind(x=c(row(mx2)), y=c(col(mx2)), z=c(mx2))
mxlr <- mxl %*%
  t(rot_z(215)) %*%
  t(rot_x(-30))

# Translate so that lowest y value is zero

mxlr[,2] <- mxlr[,2] - min(mxlr[,2])

clear3d()
spheres3d(mxlr, col=rgb(cbind(c(sh2), c(sh2), c(sh2))))

# Account for perspective; this is not correct as we're doing it purely on the
# basis of how far from the observer we are along the y axis, instead of actual
# distance.

obs <- t(matrix(
  c(
    diff(range(mxlr[,1])) / 2 + min(mxlr[,1]),
    -(diff(range(mxlr[,2])) / 2 + min(mxlr[,2])),
    diff(range(mxlr[,3])) / 2 + min(mxlr[,3])
  )
)
dist.to.obs <- sqrt(colSums((t(mxlr) - c(obs))^2))

fovh <- 60 / 180 * pi
fovv <- 60 / 180 * pi

mxlrp <- mxlr
# mxlrp[,1] <- mxlrp[,1] / ((mxlrp[,2] - obs[,2]) * tan(fovh) / 2)
# mxlrp[,3] <- mxlrp[,3] / ((mxlrp[,2] - obs[,2]) * tan(fovv) / 2)

# ggplot(as.data.frame(mxlrp)) + geom_point(aes(V1, V2, color=c(sh2)))

# Project onto a canvas, we will do this by interpreting our grid of points as a
# set of tiles where the points represent the vertices of each tile.  Start by
# transforming our point data into tile data

res.x <- 120
res.y <- 100

mxres <- cbind(mxlrp, c(sh2), dist.to.obs)
mxres[,2] <- (mxres[,2] - min(mxres[,2])) / diff(range(mxres[,2])) *
  (res.y - 1) + 1
mxres[,1] <- (mxres[,1] - min(mxres[,1])) / diff(range(mxres[,1])) *
  (res.x - 1) + 1

# Break up our data into a triangular mesh, first by doing a tile mesh, and then
# splitting each tile in two.  The tile is derived from the grid by collecting
# the top lef, top right, bottom left, and bottom right coords (these are
# represented by the 3rd dimension of the array).

nr <- nrow(mx2)
nc <- ncol(mx2)
lr.lc <- -c(seq_len(nc) * nr, seq_len(nr) + nr * (nc - 1L))
lr.fc <- -c(seq_len(nc) * nr, seq_len(nr))
fr.lc <- -c((seq_len(nc) - 1L) * nr + 1L, seq_len(nr) + nr * (nc - 1L))
fr.fc <- -c((seq_len(nc) - 1L) * nr + 1L, seq_len(nr))

tile.mesh <- c(
  mxres[lr.lc,],# drop last row, last col
  mxres[lr.fc,],# drop last row, first col
  mxres[fr.lc,],# drop first row, last col
  mxres[fr.fc,]# drop first row, first col
)
# also generate ids to track mesh elements back to original points (debugging)

mxres.id <- seq_len(nrow(mxres))
tile.mesh.id <-
  c(mxres.id[lr.lc], mxres.id[lr.fc], mxres.id[fr.lc], mxres.id[fr.fc])

dim(tile.mesh) <- c(length(tile.mesh) / (4L * ncol(mxres)), ncol(mxres), 4L)
dim(tile.mesh.id) <- c(length(tile.mesh.id) / 4L, 4L)

# now split into triangles

mesh <- array(0, c(nrow(tile.mesh)*2,ncol(tile.mesh),3))
mesh[seq_len(nrow(tile.mesh)),,] <- tile.mesh[,,1:3]
mesh[-seq_len(nrow(tile.mesh)),,] <- tile.mesh[,,2:4]
mesh.id <- matrix(0, nrow(tile.mesh.id) * 2, 3)
mesh.id[seq_len(nrow(tile.mesh.id)), ] <- tile.mesh.id[,1:3]
mesh.id[-seq_len(nrow(tile.mesh.id)), ] <- tile.mesh.id[,2:4]

# visual check that our mesh makes sense

mesh2 <- aperm(mesh, c(1, 3, 2))
mesh3 <- cbind(
  rep(seq_len(nrow(mesh2)), each=ncol(mesh2)),
  c(t(mesh2[,,1])), c(t(mesh2[,,2]))
)
# Need to reorder points in the mesh

# ggplot(as.data.frame(mesh3), aes(V2, V3, group=V1)) +
#  geom_polygon(color='grey', size=.2)

# For each triangle in the mesh, determine the set of integer x,y coordinates
# that could conceivably be in that triangle.

mesh.int <- cbind(
  id=seq_len(dim(mesh)[1]),
  xmin=pmin(mesh[,1L,1L], mesh[,1L,2L], mesh[,1L,3L]),
  xmax=pmax(mesh[,1L,1L], mesh[,1L,2L], mesh[,1L,3L]),
  ymin=pmin(mesh[,2L,1L], mesh[,2L,2L], mesh[,2L,3L]),
  ymax=pmax(mesh[,2L,1L], mesh[,2L,2L], mesh[,2L,3L])
)
# don't allow the maxs to be exactly at the boundary

xmax.bad <- mesh.int[,'xmax'] == floor(mesh.int[,'xmax'])
mesh.int[xmax.bad, 'xmax'] <- mesh.int[xmax.bad, 'xmax'] -
  (mesh.int[xmax.bad, 'xmax'] - mesh.int[xmax.bad, 'xmin']) * .01
ymax.bad <- mesh.int[,'ymax'] == floor(mesh.int[,'ymax'])
mesh.int[ymax.bad, 'ymax'] <- mesh.int[ymax.bad, 'ymax'] -
  (mesh.int[ymax.bad, 'ymax'] - mesh.int[ymax.bad, 'ymin']) * .01

mesh.int[, c('xmin', 'ymin')] <- ceiling(mesh.int[, c('xmin', 'ymin')])
mesh.int[, c('xmax', 'ymax')] <- floor(mesh.int[, c('xmax', 'ymax')])

# for each triangle, generate the set of integer coordinates that it could
# contain.

mesh.int.len.x <- pmax(mesh.int[, 'xmax'] - mesh.int[, 'xmin'] + 1, 0)
mesh.int.len.y <- pmax(mesh.int[, 'ymax'] - mesh.int[, 'ymin'] + 1, 0)
mesh.int.len <- mesh.int.len.x * mesh.int.len.y

# Use linearized coords that we'll decompose into x/y later

mesh.int.coords <- matrix(seq_len(res.x * res.y), nrow=res.x)
mesh.id.exp <- mesh.int[mesh.int.len > 0, 'id']
mesh.vals <- vector('list', length(mesh.id.exp))

for(i in seq_along(mesh.id.exp))  {
  mesh.dat <- mesh.int[mesh.id.exp[i],]
  mesh.vals[[i]] <- cbind(
    id=mesh.id.exp[i],
    coord=c(
      mesh.int.coords[
        mesh.dat['xmin']:mesh.dat['xmax'],
        mesh.dat['ymin']:mesh.dat['ymax']
] ) ) }
mesh.vals.raw <- do.call(rbind, mesh.vals)

# recover x/y integer coordinates, and compute the x for the edges

mesh.vals.id <- mesh.vals.raw[, 'id']
mesh.fin <- cbind(
  id=mesh.vals.id,
  x.int=((mesh.vals.raw[, 'coord'] - 1) %% res.x) + 1,
  y.int=((mesh.vals.raw[, 'coord'] - 1) %/% res.x) + 1
)
# compute vectors from integer points to the triangle vertices

mesh.vec <-
  mesh[mesh.vals.id,1:2,] - array(mesh.fin[,2:3], c(nrow(mesh.fin),2,3))

# compute angles between the triangle vertices

mesh.ang.1 <- acos(
  rowSums(mesh.vec[,,1] * mesh.vec[,,2]) /
  (sqrt(rowSums(mesh.vec[,,1]^2)) * sqrt(rowSums(mesh.vec[,,2]^2)))
)
mesh.ang.2 <- acos(
  rowSums(mesh.vec[,,2] * mesh.vec[,,3]) /
  (sqrt(rowSums(mesh.vec[,,2]^2)) * sqrt(rowSums(mesh.vec[,,3]^2)))
)
mesh.ang.3 <- acos(
  rowSums(mesh.vec[,,3] * mesh.vec[,,1]) /
  (sqrt(rowSums(mesh.vec[,,3]^2)) * sqrt(rowSums(mesh.vec[,,1]^2)))
)
# points that are within their triangles will have a sum of angles equal to 2pi
# Need to be a little more generous for precision pruposes, could lead to
# incorrect determinations in very corner cases

mesh.in <- mesh.fin[(mesh.ang.1 + mesh.ang.2 + mesh.ang.3 >= 2 * pi - 1e-6),]

# We need Z values and textures for each of the points we found.  Use the
# x-y distance weighted averages of the three points

mesh.dist <- sqrt(
  rowSums(
    aperm(
      (
        mesh[mesh.in[, 'id'], 1:2,] -
        array(rep(mesh.in[, 2:3], 3), c(nrow(mesh.in), 2, 3))
      ) ^ 2,
      c(1, 3, 2)
    ),
    dims=2
) )
mesh.t <- rowSums(mesh[mesh.in[, 'id'], 4,] * (1/mesh.dist)) /
  rowSums(1/mesh.dist)
mesh.d <- rowSums(mesh[mesh.in[, 'id'], 5,] * (1/mesh.dist)) /
  rowSums(1/mesh.dist)

insp <- c(6280, 6800, 9228)
insp.in <- which(mesh.in[, 'x.int'] == 88 & mesh.in[, 'y.int'] == 75)

mesh.d2 <- rowSums(mesh[insp, 5,] * (1/mesh.dist[insp.in,])) /
  rowSums(1/mesh.dist[insp.in,])

# Drop duplicates x-y values; we keep those closest to the observer; use
# untransformed

in.id <- seq_len(nrow(mesh.in))
dist.ord <- order(mesh.d)

id.unique <- -in.id[dist.ord][duplicated(mesh.in[dist.ord, 2:3])]
mesh.in.unique <- mesh.in[id.unique,]
mesh.t.unique <- mesh.t[id.unique]

mesh.in2 <- mesh.in[insp.in,]
in.id2 <- seq_len(nrow(mesh.in2))
dist.ord2 <- order(mesh.d2)
id.unique2 <- -in.id2[dist.ord2][duplicated(mesh.in2[dist.ord2, 2:3])]

id.unique <- -in.id[dist.ord][duplicated(mesh.in[dist.ord, 2:3])]
mesh.in.unique <- mesh.in[id.unique,]
mesh.t.unique <- mesh.t[id.unique]

points.in <- as.data.frame(cbind(mesh.in.unique, t=mesh.t.unique))
ggplot(as.data.frame(mesh3), aes(x=V2, y=V3, group=V1)) +
  # geom_polygon(size=.2, alpha=.4) +
  geom_point(data=points.in, aes(x=x.int, y=y.int, group=NULL, color=t))
  coord_cartesian(xlim=c(26,30), ylim=c(15,20))
#   geom_polygon(
#     data=as.data.frame(t(mesh[9718,,][1:2,])), aes(V1, V2, group=21435145123),
#     color='red', fill='blue'
#   )

  mesh[
  which(
    mesh[,1,1] <= 7 & mesh[,1,1] >= 3 &
    mesh[,2,1] <= 7 & mesh[,2,1] >= 3 &
    mesh[,1,2] <= 7 & mesh[,1,2] >= 3 &
    mesh[,2,2] <= 7 & mesh[,2,2] >= 3 &
    mesh[,1,3] <= 7 & mesh[,1,3] >= 3 &
    mesh[,2,3] <= 7 & mesh[,2,3] >= 3
  )[8]
    ,,
  ]

mesh.int.res <- cbind(id=rep(mesh.int[, 'id'], mesh.int.len), x=0, y=0)

mesh.int <- mesh.int[
  mesh.int[, 'xmin'] < mesh.int[, 'xmax'] &
  mesh.int[, 'ymin'] < mesh.int[, 'ymax'],
]

ggplot(as.data.frame(t(mesh[26,,][1:2,]))) +
geom_polygon(aes(V1, V2)

y.exp <- exp * diff(y.rng)
y.int <- seq(from=y.rng[1] - y.exp, to=y.rng[2] + y.exp, length.out=res.y)

x.rng <- range(mxlrp[,1])
x.exp <- exp * diff(x.rng)
x.int <- seq(from=x.rng[1] - x.exp, to=x.rng[2] + x.exp, length.out=res.x)

y <- findInterval(mxlrp[,2], y.int)
x <- findInterval(mxlrp[,1], x.int)

library(ggplot2)
ggplot(data.frame(x=x, y=y, z=c(sh2))) + geom_point(aes(x, y, color=z))


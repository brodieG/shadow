

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

# Account for perspective; this is not correct as we're doing it purely on the
# basis of how far from the observer we are along the y axis, instead of actual
# distance.

obs <- t(matrix(c(0, -.5 * diff(range(mxlr[,2])), mean(mxlr[,3]))))
fovh <- 60 / 180 * pi
fovv <- 60 / 180 * pi

mxlrp <- mxlr
mxlrp[,1] <- mxlrp[,1] / ((mxlrp[,2] - obs[,2]) * tan(fovh) / 2)
mxlrp[,3] <- mxlrp[,3] / ((mxlrp[,2] - obs[,2]) * tan(fovv) / 2)

ggplot(as.data.frame(mxlrp)) + geom_point(aes(V1, V2, color=c(sh2)))

# Project onto a canvas, we will do this by interpreting our grid of points as a
# set of tiles where the points represent the vertices of each tile.  Start by
# transforming our point data into tile data

res.x <- 120
res.y <- 100

mxres <- mxlrp
mxres[,2] <- (mxres[,2] - min(mxres[,2])) / diff(range(mxres[,2])) * (res.y - 1)
mxres[,1] <- (mxres[,1] - min(mxres[,1])) / diff(range(mxres[,1])) * (res.x - 1)

# Break up our data into a triangular mesh, first by doing a tile mesh, and then
# splitting each tile in two.  The tile is derived from the grid by collecting
# the top lef, top right, bottom left, and bottom right coords (these are
# represented by the 3rd dimension of the array).

nr <- nrow(mx2)
nc <- ncol(mx2)
tile.mesh <- c(
  # drop last row, last col
  mxres[-c(seq_len(nc) * nr, seq_len(nr) + nr * (nc - 1L)),],
  # drop last row, first col
  mxres[-c(seq_len(nc) * nr, seq_len(nr)),],
  # drop first row, last col
  mxres[-c((seq_len(nc) - 1L) * nr + 1L, seq_len(nr) + nr * (nc - 1L)), ],
  # drop first row, first col
  mxres[-c((seq_len(nc) - 1L) * nr + 1L, seq_len(nr)), ]
)
mesh.id <- matrix(seq_len(length(mx2)), nrow=nr)[-nr, -nc]
dim(tile.mesh) <- c(c(length(tile.mesh) / 12L), 3L, 4L)

# now split into triangles

mesh <- array(0, c(nrow(tile.mesh) * 2, 3, 3))
mesh[seq_len(nrow(tile.mesh)),,] <- tile.mesh[,,1:3]
mesh[-seq_len(nrow(tile.mesh)),,] <- tile.mesh[,,2:4]

# visual check that our mesh makes sense

mesh2 <- aperm(mesh, c(1, 3, 2))
mesh3 <- cbind(
  rep(seq_len(nrow(mesh2)), 
  each=ncol(mesh2)), c(t(mesh2[,,1])), c(t(mesh2[,,2]))
)
# Need to reorder points in the mesh

ggplot(as.data.frame(mesh3), aes(V2, V3, group=V1)) +
  geom_polygon(color='grey', size=.2)

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

mesh.int.len.x <- pmax(mesh.int[, 'xmax'] - mesh.int[, 'xmin'], 0)
mesh.int.len.y <- pmax(mesh.int[, 'ymax'] - mesh.int[, 'ymin'], 0)
mesh.int.len <- mesh.int.len.x * mesh.int.len.y

# Use linearized coords that we'll decompose into x/y later

mesh.int.coords <- matrix(seq_len(res.x * res.y), nrow=res.x)

# Special case, one point max per triangle

mesh.one <- mesh.int[mesh.int.len == 1, c('id', 'xmin', 'ymin')]

# All others

mesh.many.id <- mesh.int[mesh.int.len > 1, 'id']
mesh.many.vals <- vector('list', length(mesh.many.id))

mesh.int.res <- cbind(id=rep(mesh.int[, 'id'], mesh.int.len), x=0, y=0)

mesh.int <- mesh.int[
  mesh.int[, 'xmin'] < mesh.int[, 'xmax'] &
  mesh.int[, 'ymin'] < mesh.int[, 'ymax'],
]



y.exp <- exp * diff(y.rng)
y.int <- seq(from=y.rng[1] - y.exp, to=y.rng[2] + y.exp, length.out=res.y)

x.rng <- range(mxlrp[,1])
x.exp <- exp * diff(x.rng)
x.int <- seq(from=x.rng[1] - x.exp, to=x.rng[2] + x.exp, length.out=res.x)

y <- findInterval(mxlrp[,2], y.int)
x <- findInterval(mxlrp[,1], x.int)

library(ggplot2)
ggplot(data.frame(x=x, y=y, z=c(sh2))) + geom_point(aes(x, y, color=z))


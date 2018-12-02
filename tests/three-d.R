

# eltif <- raster::raster("~/Downloads/dem_01.tif")
# elmat1 <- matrix(
#   raster::extract(eltif,raster::extent(eltif),buffer=10000),
#   nrow=ncol(eltif),ncol=nrow(eltif)
# )
library(shadow)
# mx2 <- volcano
mx2 <- elmat1
els <- seq(-90, 90, length=25)
sun <- 180
sh2 <- ray_shade2(mx2, els, sun)
# sun <- 0
# sh2 <- rayshader::ray_shade(mx2, els, sun, lambert=FALSE)
# sh2 <- sh2[, rev(seq_len(ncol(sh2)))]
sh2 <- sh2 * .9 + .1

# Convert matrix to long format, and rotate

angle <- 3 * c(1,-1)
#rot <- rot_x(-20) %*% rot_z(215)
rot <- rot_x(-20) %*% rot_z(65)
#proj <- project_elev(mx2, sh2, rot, parallax=angle)
# system.time(
proj <- project_elev(
  mx2, sh2, rot, parallax=angle, dist=.5, resolution=c(800,800)
)
proj <- project_elev(
  mx2, sh2, diag(3), parallax=0, dist=1e6, resolution=c(800,800)
)
left <- proj[[1]]
right <- proj[[2]]


empty.rows <- which(empty_rows(left) & empty_rows(right))
empty.cols <- which(empty_cols(left) & empty_cols(right))
left <- left[-empty.rows, -empty.cols, ]
png::writePNG(left, 'persp-left.png')
right <- right[-empty.rows, -empty.cols, ]

#

left.r <- left[,,1]
left.r[left.r < 0] <- 0
left2 <- left
left2[left2 < 0] <- 0
str(as.raster(left.r))
plot(as.raster(left2))

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
res[res < 0] <- 0
png('persp-color-2.png', width=ncol(res), height=nrow(res))
par(mai=c(.25,.25,.25,.25))
par(bg='black')
plot(as.raster(res))
dev.off()

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


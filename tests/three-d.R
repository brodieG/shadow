
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

library(ggplot2)
ggplot(as.data.frame(mxlrp)) + geom_point(aes(V1, V2, color=c(sh2)))

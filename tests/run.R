## These tests are not really intended to be run on CRAN...
##
##


stop("Don't run this blindly")

eltif = raster::raster("~/Downloads/dem_01.tif")
elmat1 = matrix(
  raster::extract(eltif,raster::extent(eltif),buffer=10000),
  nrow=ncol(eltif),ncol=nrow(eltif)
)
elmat2 <- elmat1[, rev(seq_len(ncol(elmat1)))]


system.time(
yy <- rayshader::ray_shade(
  elmat2, sunangle=sun, anglebreaks=els, maxsearch=100, lambert=FALSE
)
)
system.time(
xx <- ray_shade2(elmat2, sunangle=sun, anglebreaks=els, maxsearch=100)
)
system.time(
zz <- ray_shade1(elmat1, sunangle=sun, anglebreaks=els, maxsearch=100)
)

sun <- 45
els <- seq(-90, 90, length=25)
max <- sqrt(ncol(volcano)^2+nrow(volcano)^2)

system.time(
  r0 <-
    rayshader::ray_shade(volcano, sunangle=sun, anglebreaks=els, maxsearch=max)
)
system.time(
  r1 <- ray_shade1(volcano, sunangle=sun, anglebreaks=els, maxsearch=max)
)
system.time(
  r2 <- ray_shade2(volcano, sunangle=sun, anglebreaks=els, maxsearch=max)
)

microbenchmark::microbenchmark(
  times=5,
  rayshader::ray_shade(volcano, sunangle=sun, anglebreaks=els, maxsearch=max),
  ray_shade1(volcano, sunangle=sun, anglebreaks=els, maxsearch=max),
  ray_shade2(volcano, sunangle=sun, anglebreaks=els, maxsearch=max)
)

r1 <- ray_shade1(volcano, sunangle=45, anglebreaks=els)
r2 <- ray_shade2(volcano, sunangle=45, anglebreaks=els)
system.time(r3 <- ray_shade3(elmat1, sunangle=45, anglebreaks=els))

library(ggplot2)

dims <- lapply(dim(sh.cpp), seq_len)
df <- rbind(
   cbind(do.call(expand.grid, dims), z=c(sh.cpp), type='cpp'),
   cbind(do.call(expand.grid, dims), z=c(sh.vec), type='vec'),
   cbind(do.call(expand.grid, dims), z=c(sh.for), type='for')
)

xx <- r3
xx.df <- cbind(do.call(expand.grid, lapply(dim(xx), seq_len)), z=c(xx))

ggplot(xx.df, aes(x=Var1, y=Var2, fill=z)) +
  geom_raster() 
  +
  geom_hline(yintercept=10, color='yellow') +
  geom_vline(xintercept=62, color='yellow')

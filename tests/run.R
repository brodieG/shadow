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

sun <- 45
els <- seq(-90, 90, length=25)
yy <- shadow::ray_shade2(elmat1, sunangle=sun, anglebreaks=els, maxsearch=100)

shift <- -2/180*pi
yy.r.s <- tan(shift) * elmat1 + row(elmat1)
yy.r <- apply(yy.r.s, 2, function(x) findInterval(x, seq_along(x)))
yy.r[yy.r == 0] <- 1
yy.rh <- matrix(
  yy[cbind(c(yy.r), rep(seq_len(ncol(yy)), each=nrow(yy)))], nrow(yy)
)
yy.rh.l <- yy.rh
yy.rh.r <- yy.rh

xx.r.s <- xx.r.s - floor(min(xx.r.s)) + 1L

xx.l.s <- tan(-shift) * xx + col(xx)
xx.l.s <- xx.l.s - floor(min(xx.l.s)) + 1L


xx.r <- matrix(0, nrow=nrow(xx), ncol=max(ceiling(xx.r.s)))
xx.r[

xx.l <- matrix(0, nrow=nrow(xx), ncol=max(ceiling(xx.l.s)))


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
plot_attr <- list(
  geom_raster(),
  scale_fill_gradient(low='#333333', high='#ffffff', guide=FALSE),
  ylab(NULL), xlab(NULL),
  scale_x_continuous(expand=c(0,0)),
  scale_y_continuous(expand=c(0,0)),
  theme(
    axis.text=element_text(size=6), 
    panel.spacing = unit(0, "npc"))
)
dims <- lapply(dim(yy.rh.l), seq_len)

xx.df <- rbind(
   cbind(do.call(expand.grid, dims), z=c(yy.rh.l), type='right'),
   cbind(do.call(expand.grid, dims), z=c(yy.rh.r), type='left')
)

xx <- h.c.i$z
xx.df <- cbind(do.call(expand.grid, lapply(dim(xx), seq_len)), z=c(xx))

library(ggplot2)
ggplot(xx.df, aes(x=Var1, y=Var2, fill=z)) +
  #geom_raster() + 
  plot_attr + facet_wrap(~type) +
  geom_point(
    data=data.frame(x=c(275,275), y=c(252,252),type=c('left', 'right')),
    aes(x, y, fill=NULL), shape=24, color='red'
  )


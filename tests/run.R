## These tests are not really intended to be run on CRAN...
##
##

  volcano2,

sun <- 45
els <- seq(-90, 90, length=25)
max <- sqrt(ncol(volcano)^2+nrow(volcano)^2)

ray_shade1(volcano, sunangle=sun, anglebreaks=els, maxsearch=max)
ray_shade2(volcano, sunangle=sun, anglebreaks=els, maxsearch=max)



library(ggplot2)
xx.df <- cbind(do.call(expand.grid, lapply(dim(xx), seq_len)), z=c(xx))
ggplot(xx.df, aes(x=Var1, y=Var2, fill=z)) +
  geom_raster() +
  geom_hline(yintercept=14, color='yellow') +
  geom_vline(xintercept=20, color='yellow')

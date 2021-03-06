# Shadow

## Disclaimer

This is a demo package that illustrates it is possible to write reasonably fast
ray-shading and 3D rendering code in pure R.  It will not be published to CRAN
and we have not checked whether it complies with CRAN policies or not.  The
functions included do not check their parameters and have only been lightly
tested.

## Rayshading

This package implements the ray shading functions:

* `shadow::ray_shade1`: the original for-loop based pure R shader.
* `shadow::ray_shade2`: a vectorized pure R shader
* `shadow::ray_shade3`: Like `ray_shade2`, but without bilinear interpolation

These all derive from [Tyler Morgan Wall's `rayshader`][1].  The only thing we
did is re-implement them with vectorized semantics.  `ray_shade1` is actually
taken from [Wolf Vollprecht's Next Journal article][2].


```r
elmat <- readRDS('extra/elmat.RDS') # http://tylermw.com/data/dem_01.tif.zip
system.time(shade <- shadow::ray_shade2(elmat, seq(-90, 90, length=25), 45))
```

```
##    user  system elapsed 
##   6.956   2.916  10.348
```

```r
shade.df <- cbind(do.call(expand.grid, lapply(dim(shade), seq_len)), z=c(shade))

library(ggplot2)
ggplot(shade.df, aes(x=Var1, y=Var2, fill=z)) +
  geom_raster() +
  scale_fill_gradient(low='#333333', high='#ffffff', guide=FALSE) +
  ylab(NULL) +  xlab(NULL) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text=element_text(size=6))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

## 3D Rendering

There are various functions related to 3D rendering in 'R/render.R'.  Their use
is demonstrated in the [Stereoscopy Blog
Post](https://brodieg.com/2018/12/12/three-d-pipeline/).


```r
library(shadow)
rot <- rot_x(-20) %*% rot_z(65)
rot.l <- rot %*% rot_z(2.5)
rot.r <- rot %*% rot_z(-2.5)
shadow <- ray_shade2(volcano, seq(-90, 90, length=25), sunangle=180)
elren <- mrender_elevation(
  volcano, shadow, list(rot.l, rot.r), res=1000, d=125, fov=85
)
flip <- function(x) t(x)[rev(seq_len(ncol(x))),]
elcolor <- analygraph(flip(elren[[1]]), flip(elren[[2]]))
par(bg='black', mai=numeric(4))
plot(as.raster(elcolor))
analygraph_glasses(.15)   # for the icon h/t @jurbane2
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

## Installation

This package is only available on github:

```
## devtools::install_github('brodieg/shadow'), or:
f.dl <- tempfile()
f.uz <- tempfile()
github.url <- 'https://github.com/brodieG/shadow/archive/master.zip'
download.file(github.url, f.dl)
unzip(f.dl, exdir=f.uz)
install.packages(file.path(f.uz, 'shadow-master'), repos=NULL, type='source')
unlink(c(f.dl, f.uz))
```

## Related Items

* Tyler Morgan Wall's [rayshader package][1] and [blog post][3]
* [Wolf Vollprecht's Next Journal article][2]
* My [blog post][4] on the topic.

[1]: https://github.com/tylermorganwall/rayshader
[2]: https://nextjournal.com/wolfv/how-fast-is-r-with-fastr-pythran
[3]: http://www.tylermw.com/throwing-shade/
[4]: https://www.brodieg.com/2018/10/23/do-not-shade-r/

#' Pure Base R Version of rayshader
#'
#' This is derived from the Wolf Vollprecht's
#' [adaptation](https://nextjournal.com/wolfv/how-fast-is-r-with-fastr-pythran)
#' of the original pure R implementation from Tyler Morgan Wall's [blog
#' post](http://www.tylermw.com/throwing-shade/).
#'
#' @export
#' @inheritParams rayshader::ray_shade
#' @return numeric matrix of shadow intensities between with values in
#'   &91;0,1&93;.

ray_shade2 <- function(
  heightmap, anglebreaks=seq(40, 50, 1), sunangle=315, maxsearch=100
) {
  anglebreaks <- sort(tan(anglebreaks / 180 * pi))
  sunangle <- sunangle / 180 * pi
  sinsun <- sin(sunangle)
  cossun <- cos(sunangle)

  lim.x <- if(sinsun < 0) 1 else nrow(heightmap)
  lim.y <- if(cossun < 0) 1 else ncol(heightmap)

  coords.x <- rep(seq_len(nrow(heightmap)), ncol(heightmap))
  coords.y <- rep(seq_len(ncol(heightmap)), each=nrow(heightmap))

  # For each coordinate, deterimine maximum distance towards light source before
  # we fall off grid

  dists <- as.integer(
    pmin((lim.x - coords.x) / sinsun, (lim.y - coords.y) / cossun, maxsearch)
  )
  max.dist <- max(dists)

  # We need to compute the coordinates along the path to light source.  To do
  # this generate one vector that will contain concatenated the paths for each
  # of the original points.

  sin.dists <- seq_len(max.dist) * sinsun
  cos.dists <- seq_len(max.dist) * cossun

  # coords.id tracks the origin coordinate of each element in our light path

  coords.id <- rep(seq_along(coords.x), dists)
  coords.offset <- sequence(dists)
  coords.x.id <- coords.x[coords.id]
  coords.y.id <- coords.y[coords.id]

  coords.l.x <- coords.x.id + sin.dists[coords.offset]
  coords.l.y <- coords.y.id + cos.dists[coords.offset]

  # compute angular height at each of the light coords

  heights.ang <- (
    ff_bilin(heightmap, coords.l.x, coords.l.y) -
    heightmap[cbind(coords.x.id, coords.y.id)]
  ) / coords.offset

  # compute max ang for each coord id, and then determine how many of the
  # anglebreaks are above that.

  ang.max <- tapply(heights.ang, coords.id, max)
  res <- matrix(1, nrow=nrow(heightmap), ncol=ncol(heightmap))
  res[as.integer(names(ang.max))] <- 1 - (
    findInterval(ang.max, anglebreaks, left.open=TRUE) / length(anglebreaks)
  )
  res
}
## Bilinear Interpolation
##
## A vectorized version of Wolf Vollprecht's faster_bilinear (see
## R/slow-shade.R).

ff_bilin <- function(Z, x, y) {
  i <- as.integer(x)
  j <- as.integer(y)
  XT <- x - i
  YT <- y - j

  # Expand matrix to handle OOB cases; algorithm as written uses zero
  # (implicitly), and so do we, but might be better to just use last row and
  # column values instead.  We have not tested whether the additional allocation
  # for the matrix is slower than trying to explicitly compute each of the
  # special cases.

  Z2 <- matrix(0, nrow=nrow(Z) + 1L, ncol=ncol(Z) + 1L)
  Z2[-nrow(Z2), -ncol(Z2)] <- Z

  # Create symbols for re-used vectors

  XT_1 <- 1 - XT
  YT_1 <- 1 - YT
  i1 <- i + 1L
  j1 <- j + 1L

  # compute interpolation in a vectorized manner

  (YT_1) * (XT_1) * Z2[cbind(i,  j)] +
  (YT_1) * (XT)   * Z2[cbind(i1, j)] +
  (YT)   * (XT_1) * Z2[cbind(i,  j1)] +
  (YT)   * (XT)   * Z2[cbind(i1,  j1)]
}

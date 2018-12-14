#' @export

rot_x <- function(deg) {
  r <- deg / 180 * pi
  matrix(byrow=TRUE, nrow=3, c(
         1,       0,       0,
         0,  cos(r), -sin(r),
         0,  sin(r),  cos(r))
  )
}
#' @export

rot_y <- function(deg) {
  r <- deg / 180 * pi
  matrix(byrow=TRUE, nrow=3, c(
    cos(r),       0,  sin(r),
         0,       1,       0,
   -sin(r),       0,  cos(r))
  )
}
#' @export

rot_z <- function(deg) {
  r <- deg / 180 * pi
  matrix(byrow=TRUE, nrow=3, c(
    cos(r), -sin(r),       0,
    sin(r),  cos(r),       0,
         0,       0,       1)
  )
}
#' Internally Vectorized Version of sequence
#'
#' Should be substantially faster for large `nvec`.
#'
#' @inheritParams base::sequence
#' @export

sequence2 <- function(nvec) {
  seq <- seq_len(sum(nvec))
  reset <- cumsum(nvec[-length(nvec)]) + 1L
  sub <- integer(length(seq))
  sub[reset] <- reset - 1L
  seq - cummax(sub)
}
#' @export

analygraph_glasses <- function(scale=1) {
  x.frames <- c(0, 1,  1, .6, .5, .4,  0) * .8 + .1
  y.frames <- c(1, 1, .6, .6, .7, .6, .6) * .8 + .1
  x.left <- c(.1, .4, .4, .1) * .8 + .1
  x.right <- x.left + .5 * .8
  y.left <- y.right <- c(.7, .7, .9, .9) * .8 + .1
  usr <- par()[['usr']]
  usr.x <- usr[1:2]
  usr.y <- usr[3:4]
  x <- c(x.frames, NA, x.left, NA, x.right) * diff(range(usr.x)) * scale +
    usr.x[1]
  y <- c(y.frames, NA, y.left, NA, y.right) * diff(range(usr.y)) * scale +
    (1 - scale) * diff(range(usr.y)) + usr.y[1]

  polygon(x, y, col=c('#CCCCCC', '#00FFFF', '#FF0000'), border=NA)
}
#' @export

analygraph <- function(left, right) {
  analygraph <- array(0, dim=c(dim(left), 3))
  analygraph[,,1] <- left      # red channel
  analygraph[,,2:3] <- right   # green and blue
  analygraph
}


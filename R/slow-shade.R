## A simple bilinear interpolation

faster_bilinear <- function (Z, x0, y0){
  i = floor(x0)
  j = floor(y0)
  XT = (x0 - i)
  YT = (y0 - j)
  result = (1 - YT) * (1 - XT) * Z[i, j]
  nx = nrow(Z)
  ny = ncol(Z)
  if(i + 1 <= nx){
    result = result + (1-YT) * XT * Z[i + 1, j]
  }
  if(j + 1 <= ny){
    result = result + YT * (1-XT) * Z[i, j + 1]
  }
  if(i + 1 <= nx && j + 1 <= ny){
    result = result + YT * XT * Z[i + 1, j + 1]
  }
  result
}
#' Original Pure Base R Version of rayshader
#'
#' This is based on the original Tyler Morgan Wall [blog
#' post](http://www.tylermw.com/throwing-shade/), with some clean-up to
#' functionalize, and with Wolf Vollprecht's bilinear interpolation
#' [adaptation](https://nextjournal.com/wolfv/how-fast-is-r-with-fastr-pythran)
#' function.
#'
#' inheritParams rayshader::ray_shade
#' @export
#' @return numeric matrix of shadow intensities between with values in
#'   &91;0,1&93;.

ray_shade1 <- function(
  heightmap, anglebreaks=seq(40, 50, 1), sunangle=315, maxsearch=100
) {
  elevations <- sort(tan(anglebreaks / 180 * pi))
  azimuth <- sunangle / 180 * pi
  maxdistance <- min(
    floor(sqrt(ncol(heightmap)^2 + nrow(heightmap)^2)), maxsearch
  )
  cossun <- cos(azimuth)
  sinsun <- sin(azimuth)
  heightmapshadow = matrix(1, ncol = ncol(heightmap), nrow = nrow(heightmap))

  for (i in 1:nrow(heightmap)) {
    for (j in 1:ncol(heightmap)) {
      vij = heightmap[i, j]
      for (elevation in elevations) {
        for (k in 1:maxdistance) {

          xcoord = (i + sinsun*k)
          ycoord = (j + cossun*k)

          if(xcoord > nrow(heightmap) ||
             ycoord > ncol(heightmap) ||
             xcoord < 0 || ycoord < 0) {
            break
          } else {
            tanangheight = vij + elevation * k
            if (all(c(heightmap[ceiling(xcoord), ceiling(ycoord)],
                      heightmap[floor(xcoord), ceiling(ycoord)],
                      heightmap[ceiling(xcoord), floor(ycoord)],
                      heightmap[floor(xcoord), floor(ycoord)]) < tanangheight))
              next
            if (tanangheight < faster_bilinear(heightmap, xcoord, ycoord)) {
              heightmapshadow[i, j] =  heightmapshadow[i, j] - 1 /
                length(elevations)
              break
            }
          }
        }
      }
    }
  }
  heightmapshadow
}

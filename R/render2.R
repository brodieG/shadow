# Versions of render that match exactliy to what we're doing in the blog post

#' Apply Perspective Transformation
#'
#' This function does not check its parameter.
#'
#' @export
#' @param L list of equal length numeric vectors containing named elements 'x',
#'   'y', and 'z'.  No NA elements allowed.
#' @param D scalar numeric observer distance from object along Z axis as a
#'   multiple of the object's Z-depth.
#' @return L, with 'x', and 'y' components scaled for perspective, and
#'   'x', 'y', and 'z' components normalized so that the origin is in
#'   the middle of the 'x'-'y' with the observer Z position at zero.

perspective <- function(L, D) {
  # Normalize coordinates
  L[c('x','y')] <- lapply(L[c('x','y')], function(x) x - sum(range(x)) / 2)
  z.rng <- range(L[['z']])
  z.depth <- diff(z.rng)
  L[['z']] <- L[['z']] - (z.rng[2] + D * z.depth)

  # Apply perspective transformation
  z.factor <- -1/L[['z']]
  L[c('x','y')] <- lapply(L[c('x','y')], '*', z.factor)
  L
}
#' Generate Meshes From Elevation Map
#'
#' `mesh_tile` produces a tile mesh, `mesh_tri` produces a triangle mesh.
#' The elevation map should have been converted to a list of equal length
#' vectors in "long" format.
#'
#' @rdname mesh
#' @param L list of equal length numeric vectors containing named elements 'x',
#'   'y', and 'z'.  No NA elements allowed.
#' @param dim integer(2L) the dimensions of the elevation matrix that `L`
#'   corresponds to.
#' @return a list matrix of equal length numeric vectors, with rows each
#'   representing a mesh element vertex, and as many columns as `L` has
#'   elements.

mesh_tile <- function(L, dim) {
  nr <- dim[1]
  nc <- dim[2]
  idx.raw <- matrix(seq_len(nr * nc), nr, nc)
  idx.tile <- list(
    idx.raw[-nr, -nc, drop=TRUE],
    idx.raw[-nr,  -1, drop=TRUE],
    idx.raw[ -1,  -1, drop=TRUE],
    idx.raw[ -1, -nc, drop=TRUE]
  )
  matrix(
    Map(
      function(x, y) L[[y]][idx.tile[[x]]],
      rep(seq_along(idx.tile), length(L)),
      rep(seq_along(L), each=length(idx.tile))
    ),
    nrow=length(idx.tile),
    ncol=length(L),
    dimnames=list(sprintf('v%d', seq_along(idx.tile)), names(L))
  )
}
#' @export
#' @rdname mesh

mesh_tri <- function(L, dim) {
  mesh.tile <- mesh_tile(L, dim)
  mesh.tri <- Map('c', mesh.tile[1:3,], mesh.tile[c(3,4,1),])
  dim(mesh.tri) <- c(3, ncol(mesh.tile))
  dimnames(mesh.tri) <- list(head(rownames(mesh.tile), -1), colnames(mesh.tile))
  mesh.tri
}
#' Compute Barycentric Coordinates
#'
#' @export
#' @param p a list containing equal length numeric vectors containing the
#'   coordinates of the points that we want to compute barycentric coordinates
#'   for.
#' @param v numeric list-matrix where each row corresponds to a triangle vertex
#'   data and each column the x and y coordinates for each vertex.  Each element
#'   of the underlying vectors is matched position wise to the corresponding
#'   point in `p`, so all vectors in `p` and `v` are expected to be the same
#'   length.
#' @return a list with there numeric vectors containing the barycentric
#'   coordinates.

barycentric <- function(p, v) {
  det <-  (v[[2,'y']]-v[[3,'y']])*(v[[1,'x']]-v[[3,'x']]) +
          (v[[3,'x']]-v[[2,'x']])*(v[[1,'y']]-v[[3,'y']])

  l1 <- (
          (v[[2,'y']]-v[[3,'y']])*(  p[['x']]-v[[3,'x']]) +
          (v[[3,'x']]-v[[2,'x']])*(  p[['y']]-v[[3,'y']])
        ) / det
  l2 <- (
          (v[[3,'y']]-v[[1,'y']])*(  p[['x']]-v[[3,'x']]) +
          (v[[1,'x']]-v[[3,'x']])*(  p[['y']]-v[[3,'y']])
        ) / det
  l3 <- 1 - l1 - l2
  list(l1, l2, l3)
}

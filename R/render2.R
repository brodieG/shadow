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

mesh_tri <- function(L, dim, order=FALSE) {
  mesh.tile <- mesh_tile(L, dim)
  mesh.tri <- Map('c', mesh.tile[1:3,], mesh.tile[c(3,4,1),])
  dim(mesh.tri) <- c(3, ncol(mesh.tile))
  dimnames(mesh.tri) <- list(head(rownames(mesh.tile), -1), colnames(mesh.tile))

  if(order) {
    zord <- order(Reduce('+', mesh.tri[,'z']))
    mesh.tri[['z']] <- NULL   # don't need z val anymore
    mesh.tri[] <- lapply(mesh.tri, '[', zord)
  }
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

to_long <- function(x) rbind(x=c(row(x)), y=c(col(x)), z=c(x))
rotate <- function(x, rot) rot %*% x
mx_as_list <- function(x, texture)
  c(
    setNames(lapply(seq_len(nrow(x)), function(i) x[i,]), c('x','y','z')),
    t=list(c(texture))
  )

perspective <- function(x, d) {
  z.rng <- range(x[['z']])
  ZD <- diff(z.rng)

  ## Normalize XY coordinates so observer is centered in midpoint of x-y extent
  ## and at Z=0
  xp <- x
  xp[c('x','y')] <- lapply(x[c('x','y')], function(x) x - sum(range(x)) / 2)
  xp[['z']] <- xp[['z']] - (z.rng[2] + d * ZD)

  ## Apply perspective transformation
  z.factor <- -1 / xp[['z']]
  xp[c('x','y')] <- lapply(xp[c('x','y')], '*', z.factor)
  xp
}
bounding_boxes <- function(mesh) {
  mins <- apply(
    mesh[,c('x','y')], 2,
    function(x) list(as.integer(ceiling(do.call(pmin, x))))
  )
  maxs <- apply(
    mesh[,c('x','y')], 2,
    function(x) list(as.integer(do.call(pmax, x)))
  )
  list(mins=lapply(mins, '[[', 1L), maxs=lapply(maxs, '[[', 1L))
}
candidate_pixels <- function(bb) {
  mins <- bb[['mins']]
  maxs <- bb[['maxs']]
  dims.x <- pmax(maxs[['x']] - mins[['x']] + 1, 0)
  dims.y <- pmax(maxs[['y']] - mins[['y']] + 1, 0)
  dims.xy <- dims.x * dims.y
  bb.id <- rep(seq_along(dims.x), dims.xy) # if dims.xy zero, then id omitted

  dims.val <- dims.xy > 0
  dims.x.val <- dims.x[dims.val]
  dims.y.val <- dims.y[dims.val]
  x.lens <- rep(dims.x.val, dims.y.val)
  x.off <- sequence2(x.lens) - 1L
  y.off <- rep(sequence2(dims.y.val), x.lens) - 1L

  ## Add back the min values of the bounding boxes
  p.x <- x.off + mins[['x']][bb.id]
  p.y <- y.off + mins[['y']][bb.id]
  list(id=bb.id, x=p.x, y=p.y)
}
mesh_expand <- function(mesh, ids) {
  mesh.all <- mesh
  mesh.all[] <- lapply(mesh, '[', ids)
  mesh.all
}
shade <- function(texture, bc) Reduce('+', Map('*', texture, bc))

rasterize <- function(mesh, resolution) {
  bb <- bounding_boxes(mesh)
  p.cand <- candidate_pixels(bb)
  mesh.all <- mesh_expand(mesh, p.cand[['id']])
  bc <- barycentric(p.cand, mesh.all)
  inbounds <- Reduce('&', lapply(bc, '>=', 0))

  texture <-
    shade(lapply(mesh.all[,'t'], '[', inbounds), lapply(bc, '[', inbounds))

  p.in <- lapply(p.cand[c('x','y')], '[', inbounds)
  p.raster <- matrix(numeric(0), resolution[1], resolution[2])
  p.raster[do.call(cbind, unname(p.in))] <- texture
  p.raster
}

project_and_scale <- function(elevation, texture, rotation, resolution, d) {
  el.long <- to_long(elevation)
  r <- rotate(el.long, rotation)
  rl <- mx_as_list(r, texture)
  rlp <- perspective(rl, d)

  x.rng <- range(rlp[['x']])
  y.rng <- range(rlp[['y']])
  asp <- diff(x.rng) / diff(y.rng)
  y.res <- as.integer(round(resolution * asp))

  rlp[['x']] <- (rlp[['x']] - x.rng[1]) / diff(x.rng) * (resolution - 1) + 1
  rlp[['y']] <- (rlp[['y']] - y.rng[1]) / diff(y.rng) * (y.res - 1) + 1
  attr(rlp, 'resolution') <- as.integer(c(x=resolution, y=y.res))
  rlp
}
#' Convert Elevation Map as Triangle Polygon Mesh
#'
#' @export

elevation_as_mesh <- function(elevation, texture, rotation, resolution, d) {
  rlp <- project_and_scale(elevation, texture, rotation, resolution, d)
  mesh_tri(rlp, dim(elevation), order=TRUE)
}

#' Render a 3D Elevation as a 2D Image
#'
#' @param elevation numeric matrix of elevations; each matrix element is assumed
#'   to be equally spaced in both X and Y dimensions.
#' @param texture a numeric matrix with values in [0,1] of the same dimensions
#'   as `elevation`, where 1 means light and 0 means dark.  Currently only
#'   grayscale textures are supported.
#' @param 3 x 3 numeric rotation matrix.
#' @param resolution integer(1L) the width of the output in pixels assuming
#'   there is no rotation.  Actual width will depend on whether rotation causes
#'   the width to change.
#' @param d numeric(1L) distance of the observer from nearest part of the model
#'   as a multiple of the total depth of the model along the observation axis.
#' @export

render_elevation <- function(elevation, texture, rotation, resolution, d) {
  stopifnot(identical(dim(elevation), dim(texture)))
  resolution <- as.integer(resolution)
  rlp <- project_and_scale(elevation, texture, rotation, resolution, d)
  mesh <- mesh_tri(rlp, dim(elevation), order=TRUE)
  rasterize(mesh, attr(rlp, 'resolution'))
}

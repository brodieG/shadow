# Versions of render that match exactliy to what we're doing in the blog post

#' Apply Relative Perspective Transformation
#'
#' Observer position is specified as a multiple of the Z depth of the object,
#' measured from the nearest Z point of the object.
#'
#' @export
#' @param L list of equal length numeric vectors containing named elements 'x',
#'   'y', and 'z'.  No NA elements allowed.
#' @param d scalar numeric observer distance from object along Z axis as a
#'   multiple of the object's Z-depth, or for `persp_abs` the distance from the
#'   midpoint of the model..
#' @param mode character in "rel", "abs". If "rel" then `d` is interpreted as
#'   multiples of Z depth.  If "abs" it is taken directly as the distance from
#'   the zero Z value, which should be the point around which the model was
#'   rotated.
#' @return L, with 'x', and 'y' components scaled for perspective, and
#'   'x', 'y', and 'z' components normalized so that the origin is in
#'   the middle of the 'x'-'y' with the observer Z position at zero.

persp_rel <- function(L, d) {
  ## Normalize XY coordinates so observer is centered in midpoint of x-y extent
  ## and at Z=0
  LP <- L
  LP[c('x','y','z')] <-
      lapply(LP[c('x','y','z')], function(x) x - sum(range(x)) / 2)
  z.rng <- range(LP[['z']])
  ZD <- diff(z.rng)
  LP[['z']] <- LP[['z']] - (z.rng[2] + d * ZD)

  ## Apply perspective transformation
  z.factor <- -1 / LP[['z']]
  LP[c('x','y')] <- lapply(LP[c('x','y')], '*', z.factor)
  LP
}
#' @export
#' @rdname persp_rel

persp_abs <- function(L, d, fov) {
  L[['z']] <- L[['z']] - d
  fovv <- tan(fov / 2 / 180 * pi)
  L[c('x','y')] <- lapply(L[c('x','y')], function(y) y / (fovv * -L[['z']]))
  L
}
#'
#' Generate Meshes From Elevation Map
#'
#' `mesh_tile` produces a tile mesh, `mesh_tri` produces a triangle mesh.
#' The elevation map should have been converted to a list of equal length
#' vectors in "long" format.
#'
#' @export
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

to_long <- function(elevation) {
  # rbind(x=c(row(elevation)), y=c(col(elevation)), z=c(elevation))
  x <- c(row(elevation))
  y <- c(col(elevation))
  rbind(
    x=x - sum(range(x)) / 2,
    y=y - sum(range(y)) / 2,
    z=c(elevation - sum(range(elevation)) / 2)
  )
}
mx_as_list <- function(x, texture)
  c(
    setNames(lapply(seq_len(nrow(x)), function(i) x[i,]), c('x','y','z')),
    t=list(c(texture))
  )

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

# Rasterize Triangle Mesh

rasterize <- function(mesh, resolution, zord, empty) {
  bb <- bounding_boxes(mesh)
  p.cand <- candidate_pixels(bb)
  mesh.all <- mesh_expand(mesh, p.cand[['id']])
  bc <- barycentric(p.cand, mesh.all)
  inbounds <- Reduce('&', lapply(bc, function(x) !is.na(x) & x >= 0))
  bc.in <- lapply(bc, '[', inbounds)
  texture <- shade(lapply(mesh.all[,'t'], '[', inbounds), bc.in)
  p.in <- lapply(p.cand[c('x','y')], '[', inbounds)

  if(zord == 'pixel') {
    zord  <- order(shade(lapply(mesh.all[,'z'], '[', inbounds), bc.in))
    p.in <- lapply(p.in, '[', zord)
    texture <- texture[zord]
  }

  p.raster <- matrix(empty, resolution[1], resolution[2])
  p.raster[do.call(cbind, unname(p.in))] <- texture
  p.raster
}
rotate <- function(elevation, texture, rotation) {
  el.long <- to_long(elevation)
  r <- rotation %*% el.long
  rl <- mx_as_list(r, texture)
}
# Scale to Pixel Coordinates
#
# `scale_abs` is used when perspective was applied on an absolute basis with a
# field of view and may lead to either cropping (well, bunching really) or empty
# space.  You may need to change the `fov` value to get a good crop.
#
# `scale_rel` will just fit to a particular pixel width.

scale_abs <- function(L, resolution) {
  x.rng <- range(L[['x']])
  y.rng <- range(L[['y']])

  scale_res <- function(w, res) {
    w <- w * res + res / 2
    w[w < 1] <- 1
    w[w > res] <- res
    w
  }
  L[c('x','y')] <- Map(scale_res, L[c('x','y')], c(resolution, resolution))
  L
}
scale_rel <- function(L, resolution) {
  x.rng <- range(L[['x']])
  y.rng <- range(L[['y']])
  asp <- diff(x.rng) / diff(y.rng)
  y.res <- as.integer(round(resolution / asp))

  L[['x']] <- (L[['x']] - x.rng[1]) / diff(x.rng) * (resolution - 1) + 1
  L[['y']] <- (L[['y']] - y.rng[1]) / diff(y.rng) * (y.res - 1) + 1
  attr(L, 'resolution') <- as.integer(c(x=resolution, y=y.res))
  L
}
#' Convert Elevation Map as Triangle Polygon Mesh
#'
#' @export

elevation_as_mesh <- function(elevation, texture, rotation, d) {
  rl <- rotate(elevation=elevation, texture=texture, rotation=rotation)
  rlp <- persp_rel(rl, d)
  mesh <- mesh_tri(rlp, dim(elevation), order=TRUE)
}

#' Render a 3D Elevation as a 2D Image
#'
#' `mrender_elevation` runs multiple projections with the same perspective and
#' field of view parameters so that they can be used together in a video.
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

render_elevation <- function(
  elevation, texture, rotation, resolution, d=1, zord='mesh', empty=0
) {
  stopifnot(
    identical(dim(elevation), dim(texture)),
    isTRUE(zord %in% c('mesh', 'pixel'))
  )
  rl <- rotate(elevation=elevation, texture=texture, rotation=rotation)
  rlp <- persp_rel(rl, d)
  resolution <- as.integer(resolution)
  rlps <- scale_rel(rlp, resolution)
  mesh <- mesh_tri(rlps, dim(elevation), order=zord == 'mesh')
  rasterize(mesh, attr(rlps, 'resolution'), zord, empty)
}

#' @export

mrender_elevation <- function(
  elevation, texture, rotations, resolution, d=diff(range(elevation)) * 2,
  fov=60, zord='mesh', empty=0, trim=TRUE
) {
  if(!is.list(elevation)) elevation <- list(elevation)
  if(!is.list(texture)) texture <- list(texture)
  if(!is.list(d)) d <- list(d)
  stopifnot(is.list(rotations))

  elevation <- rep(elevation, length.out=length(rotations))
  texture <- rep(texture, length.out=length(rotations))
  d <- rep(d, length.out=length(rotations))

  stopifnot(
    identical(lapply(elevation, dim), lapply(texture, dim)),
    tail(duplicated(t(sapply(elevation, dim))), -1L),
    isTRUE(zord %in% c('mesh', 'pixel'))
  )
  #Perspective transform for each resolution

  resolution <- as.integer(resolution)
  RL <- Map(rotate, elevation=elevation, texture=texture, rotations)
  RLP <- Map(persp_abs, RL, fov=fov, d=d)

  # Need to adjust resolutions so that they mean the same thing across the
  # various frames.

  RLPS <- Map(scale_abs, RLP, resolution)
  MESH <- lapply(RLPS, mesh_tri, dim(elevation[[1]]), order=zord == 'mesh')
  ELREN <- mapply(
    rasterize, MESH,
    MoreArgs=list(zord=zord, empty=empty, resolution=c(resolution, resolution)),
    SIMPLIFY=FALSE
  )
  # Remove totally blank lines; this currently only works if `empty` is set to 0
  if(trim) {
    cols <- rle(Reduce('+', lapply(ELREN, colSums)))
    rows <- rle(Reduce('+', lapply(ELREN, rowSums)))
    col.drop <- c(
      if(!cols[['values']][1]) -seq_len(cols[['lengths']][1]),
      if(!cols[['values']][length(cols[['values']])])
        -seq(
          from=ncol(ELREN[[1]]),
          length.out=cols[['lengths']][length(cols[['values']])],
          by=-1
        )
    )
    row.drop <- c(
      if(!rows[['values']][1]) -seq_len(rows[['lengths']][1]),
      if(!rows[['values']][length(rows[['values']])])
        -seq(
          from=ncol(ELREN[[1]]),
          length.out=rows[['lengths']][length(rows[['values']])],
          by=-1
        )
    )
    if(!length(col.drop)) col.drop <- seq_len(ncol(ELREN[[1]]))
    if(!length(row.drop)) row.drop <- seq_len(nrow(ELREN[[1]]))
    lapply(ELREN, '[', row.drop, col.drop)
  } else ELREN
}

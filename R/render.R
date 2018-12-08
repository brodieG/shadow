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
#' Create Triangular Mesh From Elevation Data
#'
#' Elevation data should already be in long format, and ordered by row and
#' column.
#'
#' @export
#' @param dat numeric matrix with four columns representing respectively the x,
#'   y, z, and texture values of the elevation map.
#' @param nr number of rows in original elevation matrix
#' @param nc number of columns in original elevation matrix
#' @return 3 x 5 list matrix where rows index the vertex, and the columns
#'   represent the x, y, z, and texture values.

triangular_mesh <- function(dat, nr, nc) {
  stopifnot(identical(dim(dat)[2], 4L), nr * nc == nrow(dat))

  # create a tile mesh by combining four copies of our original matrix,
  # dropping first/last row/col in turn.  The arithmetic is to map location
  # of first/last row/col to our long input matrix `dat`.

  lr.lc <- -c(seq_len(nc) * nr, seq_len(nr) + nr * (nc - 1L))
  lr.fc <- -c(seq_len(nc) * nr, seq_len(nr))
  fr.lc <- -c((seq_len(nc) - 1L) * nr + 1L, seq_len(nr) + nr * (nc - 1L))
  fr.fc <- -c((seq_len(nc) - 1L) * nr + 1L, seq_len(nr))

  tile.mesh <- c(
    dat[lr.lc,],   # drop last row, last col
    dat[lr.fc,],   # drop last row, first col
    dat[fr.lc,],   # drop first row, last col
    dat[fr.fc,]    # drop first row, first col
  )
  # store mesh data in a 3D array, where the indices of the 3rd dimension
  # represent which vertex of the mesh element is described.

  dim(tile.mesh) <- c(length(tile.mesh) / (4L * ncol(dat)), ncol(dat), 4L)

  # now split tiles into triangles

  mesh <- array(0, c(nrow(tile.mesh)*2,ncol(tile.mesh),3))
  mesh[seq_len(nrow(tile.mesh)),,] <- tile.mesh[,,1:3]
  mesh[-seq_len(nrow(tile.mesh)),,] <- tile.mesh[,,2:4]
  dimnames(mesh) <- list(NULL, c('x', 'y', 'z', 't'), NULL)

  # convert to a list matrix, where the rows index the vertices and the columns
  # the data.  We do this because we will be subsetting by columns a lot

  mesh.list.dim <- rev(dim(mesh)[-1])
  list.ids <- do.call(expand.grid, lapply(mesh.list.dim, seq_len))
  mesh.list <- Map(
    function(vertex, col) mesh[,col,vertex,drop=TRUE],
    list.ids[[1]], list.ids[[2]]
  )
  dim(mesh.list) <- mesh.list.dim
  dimnames(mesh.list) <- rev(dimnames(mesh)[-1])
  mesh.list
}
#' Generate Candidate Grid Points Within Triangles
#'
#' For each triangle in the mesh, determine which points in the integer grid
#' could potentially be in the mesh
#'
#' @export
#' @param list-matrix of numeric vectors

candidate_points <- function(mesh) {
  # assume all points between max and min on y/x could be in mesh
  points.int <- list(
    id=seq_along(mesh[[1,'x']]),
    xmin=do.call(pmin, mesh[,'x']), xmax=do.call(pmax, mesh[,'x']),
    ymin=do.call(pmin, mesh[,'y']), ymax=do.call(pmax, mesh[,'y'])
  )
  # don't allow the maxs to be exactly at the boundary

  xmax.bad <- points.int[['xmax']] == floor(points.int[['xmax']])
  points.int[['xmax']][xmax.bad] <- with(
    lapply(points.int, '[', xmax.bad), xmax - (xmax  - xmin) * .01
  )
  ymax.bad <- points.int[['ymax']] == floor(points.int[['ymax']])
  points.int[['ymax']][ymax.bad] <- with(
    lapply(points.int, '[', ymax.bad), ymax - (ymax  - ymin) * .01
  )
  xymin <- c('xmin', 'ymin')
  xymax <- c('xmax', 'ymax')
  points.int[xymin] <-
    lapply(points.int[xymin], function(x) as.integer(ceiling(x)))
  points.int[xymax] <- lapply(points.int[xymax], as.integer)

  # for each triangle, generate the set of integer coordinates that it could
  # contain.

  p.int.len.x <- pmax(points.int[['xmax']] - points.int[['xmin']] + 1L, 0L)
  p.int.len.y <- pmax(points.int[['ymax']] - points.int[['ymin']] + 1L, 0L)
  p.int.len <- p.int.len.x * p.int.len.y

  # For each id, need the equivalent of `expand.grid` of the xmin:xmax and
  # ymin:ymax sequences.  The following logic replicates that by doing it for
  # 1:(xmax-xmin) and 1:(ymax-ymin); we add back xmin and xmax later.  Contorted
  # logic is to minimize calls to R functions as this is bottleneck.

  seq.id <- rep(points.int[['id']], p.int.len)
  y.lens <- rep(p.int.len.y, p.int.len.x)

  seq.x <- rep(sequence(p.int.len.x), y.lens) - 1L
  seq.x.delta <- cumsum(y.lens[-length(y.lens)]) + 1L

  seq.y <- integer(length(seq.x))
  seq.y[seq.x.delta] <- seq.x.delta - 1L
  seq.y <- seq_along(seq.y) - cummax(seq.y) - 1L

  # add back xmin and xmax

  list(
    id=seq.id,
    x=seq.x + points.int[['xmin']][seq.id],
    y=seq.y + points.int[['ymin']][seq.id]
  )
}
#' Compute Barycentric Coordinates
#'
#' Taken from [wikipedia](https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Conversion_between_barycentric_and_Cartesian_coordinates).
#'
#' @export
#' @note all the input vectors are expected to be of the same length.
#' @param x,y, numeric vectors with cartesian coordinates of the points to
#'   compute barycentric coordinates for.
#' @param x1, x2, x3, y1, y2, y3, numeric vectors of the cartesian
#'   coordinates of the vertices of the triangles.
#'   to matrix with two columns 'x' and 'y':
#'   the points to compute
#'   the barycentric coordinates for.

barycentric_coord <- function(x, x1, x2, x3, y, y1, y2, y3) {
  y2.y3 <- y2 - y3
  x3.x2 <- x3 - x2
  x1.x3 <- x1 - x3
  x.x3 <- x - x3
  y.y3 <- y - y3
  det <- y2.y3 * x1.x3 + x3.x2 * (y1 - y3)

  l1 <- (y2.y3 * x.x3 + x3.x2 * y.y3) / det
  l2 <- ((y3 - y1) * x.x3 + x1.x3 * y.y3) / det
  l3 <- 1 - l1 - l2

  list(l1, l2, l3)
}
#' Compute Texture of Each Pixels
#'
#' Use barycentric coordinates to compute the vertex weighted texture of each
#' pixel within each mesh triangle.  We also compute the z value of each pixel
#' in the same manner.  Finally, we eliminate any pixels not in the mesh
#' triangle.
#'
#' @export
#' @param p numeric matrix with 'x' and 'y' coords, and an 'id' of what mesh
#'   triangle the point corresponds to in theory.
#' @param m numeric array with 'x', 'y', 'z', and 't' values for each
#'   vertex.

points_texture <- function(p, m) {
  # expand mesh data to line up with our points
  M <- lapply(m, '[', p[['id']])
  dim(M) <- dim(m)
  dimnames(M) <- dimnames(m)

  # compute barycentric weighted texture and z values
  bary <- do.call(
    barycentric_coord, unname(c(p['x'], M[,'x'], p['y'], M[,'y']))
  )

  # drop pixes outside of mesh triangles.
  inbounds <- bary[[1]] >= 0 & bary[[2]] >= 0 & bary[[3]] >= 0
  bary.in <- lapply(bary, '[', inbounds)
  M.in <- lapply(M[,c('t', 'z')], '[', inbounds)
  p.in <- lapply(p, '[', inbounds)
  dim(M.in) <- c(3L, 2L)
  dimnames(M.in) <- list(NULL, c('t', 'z'))

  # compute point texture and z vals
  p.t <- Reduce('+', Map('*', bary.in, M.in[,'t']))
  p.z <- Reduce('+', Map('*', bary.in, M.in[,'z']))
  c(p.in, list(z=p.z, t=p.t))
}
#' @export

empty_rows <- function(arr) (rowSums(arr) == -3 * ncol(arr))

#' @export

empty_cols <- function(arr) (rowSums(colSums(arr)) == -3 * nrow(arr))

#' Projection
#'
#' @param resolution integer(2L) output resolution where the first value is the
#'   x resolution, and the second the y resolution.
#' @param parallax numeric vector of degrees to add to the z angle of
#'   `rotation`.  A view will be rendered for each angle and returned in the
#'   result list.  This makes it easy to generate stereoscopic images.
#' @param dist numeric(1L) distance of the observer from the x-y plane as a
#'   multiple of the z distance spanned by the elevation surface after rotation.
#'   currently the observer is always centered on the x-y plane.
#' @param fov numeric(1L) field of view in degrees.
#' @export

project_elev <- function(
  elevation, texture, rotation, parallax=0, dist=2,
  resolution=as.integer(dim(elevation) * 1.2)
) {
  mxl <- rbind(c(row(elevation)), c(col(elevation)), c(elevation))
  res <- vector('list', length(parallax))

  for(i in seq_along(parallax)) {
    mxlr <- rotation %*% rot_z(parallax[i]) %*% (mxl)
    ranges <- apply(mxlr, 1, range)
    obs <- c(
      diff(ranges[,1])/2 + ranges[1,1],
      diff(ranges[,2])/2 + ranges[1,2],
      dist * diff(ranges[,3]) + ranges[2,3]
    )
    mxlro <- mxlr - obs
    dimnames(mxlro)[[1]] <- c('x', 'y', 'z')

    # apply parallax rotation and compute observer position
    # compute projection

    mxlrpp <- mxlro
    mxlrpp['x',] <- -1/mxlro['z',] * mxlro['x',]
    mxlrpp['y',] <- -1/mxlro['z',] * mxlro['y',]

    # df <- as.data.frame(cbind(t(mxlrpp), color=c(sh2)))
    # print(ggplot(df[order(df$V3),]) + geom_point(aes(V1, V2, color=color)))
    # print(ggplot(df[order(df$V3),]) + geom_point(aes(V1, V2, color=V3)))
    # mxlrp[,3] <- mxlrp[,3] / ((mxlrp[,2] - obs[,2]) * tan(fovv) / 2)

    # ggplot(as.data.frame(t(mxlr))) + geom_point(aes(V1, V2, color=c(sh2)))
    # ggplot(as.data.frame(t(mxlrpp))) + geom_point(aes(V1, V2, color=c(sh2)))
    # ggplot(as.data.frame(cbind(t(mxl[1:2,]), z=c(texture)))) +
    #   geom_point(aes(x, y, color=z))

    # Turn to list format, rescale, and add meta data

    lres <- c(
      lapply(seq_len(nrow(mxlrpp)), function(i) mxlrpp[i,]), list(texture)
    )
    names(lres) <- c('x', 'y', 'z', 't')

    lres[['y']] <-
      (lres[['y']] - min(lres[['y']])) / diff(range(lres[['y']])) *
      (resolution[2] - 1) + 1
    lres[['x']] <-
      (lres[['x']] - min(lres[['x']])) / diff(range(lres[['x']])) *
      (resolution[1] - 1) + 1

    # print(ggplot(as.data.frame(mxres)) + geom_point(aes(V1, V2, color=texture)))

    # Generate coordinates for triangular mesh from our points.

    mesh <- mesh_tri(lres, dim(elevation))

    # mesh2 <- aperm(mesh, c(1, 3, 2))
    # mesh3 <- cbind(
    #   rep(seq_len(nrow(mesh2)), each=ncol(mesh2)),
    #   c(t(mesh2[,,1])), c(t(mesh2[,,2])), c(t(mesh2[,,5]))
    # )
    # meshdf <- as.data.frame(mesh3)
    # meshdf$z <- ave(meshdf$V4, meshdf$V1, FUN=mean)
    # meshdf$id <- rep(seq_len(nrow(meshdf)/3), each=3)
    # meshdf$id <- factor(meshdf$id, levels=sort(unique(meshdf$id), dec=T))
    # ggplot(meshdf, aes(V2, V3, group=id)) +
    #  geom_polygon(color='grey', size=.05)

    points.raw <- candidate_points(mesh)
    points.dat <- points_texture(points.raw, mesh)
    # resolve duplicates by putting nearest points last
    points.ord <- order(points.dat[['z']])
    points <- lapply(points.dat, '[', points.ord)

    # ggplot(as.data.frame(points.dat)) +
    # geom_point(aes(x=x.int, y=y.int, color=texture))

    res.mx <- matrix(-1, nrow=resolution[1], ncol=resolution[2])
    res.mx[do.call(cbind, points[c('x', 'y')])] <- points[['t']]

    # transformations so renders okay as raster / png as those plot the matrix
    # in the same way it displays in the terminal (rows vertically which and in
    # increasing index value), as opposed to in traditional x y plot like RGL.

    res.mx <- t(res.mx)
    res.mx <- res.mx[rev(seq_len(nrow(res.mx))), ]
    res[[i]] <- array(res.mx, c(dim(res.mx), 3L))
  }
  res
}


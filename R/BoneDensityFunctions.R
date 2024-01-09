#library(tidyverse)
#library(rdist)
#library(oro.nifti)
#library(rgl)
#library(Rvcg)
#library(ptinpoly)


#' check landmarks are close to the bone
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh_path String. Filepath to triangulated surface mesh in
#' ply or stl format
#' @param landmark_path String. Filepath to landmark data in .fcsv format from
#' 3D Slicer
#' @param threshold Numeric. Distance landmark can be from surface without
#' warning being thrown
#' @return String. Returns a message warning that landmarks are not on bone
#' surface
#' @importFrom Rvcg vcgImport
#' @importFrom rdist cdist
#' @importFrom utils read.csv
#' @export
landmark_check <- function(surface_mesh_path, landmark_path, threshold = 1.0) {
  surface_mesh <- vcgImport(surface_mesh_path)
  vertices <- t(surface_mesh$vb)[, c(1:3)]
  landmarks <- read.csv(landmark_path, skip = 3, header = FALSE)[2:4]

  dists <- c()
  for (i in 1:nrow(landmarks)) {
    x <- cdist(vertices, landmarks[i, ])
    dists <- c(dists, min(x))
  }

  # return message if landmarks not on bone surface
  if (any(dists > threshold)) {print("landmarks not on bone surface")}
}


#' Sigma beta CT calculations
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param table_height Numeric
#' @param calibration_curves matrix
#' @param scanner String "CT1" or "CT2"
#' @param return_coeff String
#' @return Vector with two elements, sigma and beta
#' @importFrom dplyr filter
#' @importFrom stats approx
#' @importFrom magrittr "%>%"
#' @export
ct_coefficients <- function(table_height, calibration_curves, scanner, return_coeff = "sigma") {
  Scanner <- NULL

  # get curves for scanner
  calibration_curves <- calibration_curves %>% filter(Scanner == scanner)
  sigmaCT <- approx(calibration_curves$TableHeight, calibration_curves$sigma, table_height)$y
  betaCT <- approx(calibration_curves$TableHeight, calibration_curves$beta, table_height)$y

  # return
  if(return_coeff == "sigma") {return(sigmaCT)}
  if(return_coeff == "beta") {return(betaCT)}
}


#' Redefine surface points
#' @author Scott Telfer \email{scott.telfer@gmail.com} Adapted from geomorph
#' @param surface_mesh Mesh file
#' @param landmarks Data frame. Contains 3D coords of landmarks
#' @param no_surface_sliders Numeric. Number of surface points to add
#' @return Data frame. 3D coords of remapped surface points
#' @importFrom stats kmeans
#' @export
surface_points_template <- function(surface_mesh, landmarks, no_surface_sliders) {
  # get points from surface
  if (is.list(surface_mesh)) {
    vertices <- t(surface_mesh$vb)[, c(1:3)]
  } else {
    vertices <- surface_mesh
  }
  colnames(vertices) <- c("xpts", "ypts", "zpts")

  # landmarks
  lmk.add <- NULL
  for(i in 1:nrow(landmarks)){
    lmk.add <- rbind(lmk.add, which.min(sqrt((landmarks[i, 1] - vertices[, 1]) ^ 2 +
                                             (landmarks[i, 2] - vertices[, 2]) ^ 2 +
                                             (landmarks[i, 3] - vertices[, 3]) ^ 2))[1])}
  nlandmarks <- nrow(landmarks)
  vertices <- vertices[-lmk.add, ]

  # calculate new points
  colnames(landmarks) <- c("xpts", "ypts", "zpts")
  new_surface_points <- rbind(landmarks,
                              kmeans(x = vertices, centers = no_surface_sliders,
                                     iter.max = 25)$centers)

  # return
  return(new_surface_points)
}


#' New surface points from template
#' @author Scott Telfer \email{scott.telfer@gmail.com} Adpated from geomorph
#' @param surface_mesh List. Mesh data imported via ply_import function
#' @param landmarks Data frame. Contains 3D coords of landmarks
#' @param template Data frame. 3D coords of remapped surface points
#' @return Data frame. 3D coords of remapped surface points
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats dist
#' @export
surface_points_new <- function(surface_mesh, landmarks, template) {
  ## helper functions
  rotate.mat <- function(M, Y){
    k <- ncol(M)
    M <- cs.scale(M); Y <- cs.scale(Y)
    MY <- crossprod(M, Y)
    sv <- La.svd(MY, k, k)
    u <- sv$u; u[, k] <- u[, k] * determinant(MY)$sign
    v <- t(sv$vt)
    tcrossprod(v, u)
  }

  csize <- function(x) sqrt(sum(center(as.matrix(x))^2))

  center <- function(x){
    if(is.vector(x)) x - mean(x) else {
      x <- as.matrix(x)
      dims <- dim(x)
      fast.center(x, dims[1], dims[2])
    }
  }

  fast.center <- function(x, n, p){
    m <- colMeans(x)
    x - rep.int(m, rep_len(n, p))
  }

  cs.scale <- function(x) x/csize(x)

  fast.solve <- function(x) {
    x <- as.matrix(x)
    if(det(x) > 1e-8) {
      res <- try(chol2inv(chol(x)), silent = TRUE)
      if(class(res) == "try-error") res <- fast.ginv(x)
    } else res <- fast.ginv(x)
    return(res)
  }

  tps2d3d <- function(M, matr, matt, PB = TRUE){		#DCA: altered from J. Claude 2008
    p <- dim(matr)[1]; k <- dim(matr)[2]; q <- dim(M)[1]
    Pdist <- as.matrix(stats::dist(matr))
    ifelse(k == 2, P <- Pdist^2*log(Pdist^2), P <- Pdist)
    P[which(is.na(P))] <- 0
    Q <- cbind(1, matr)
    L <- rbind(cbind(P, Q), cbind(t(Q), matrix(0, k + 1, k + 1)))
    m2 <- rbind(matt, matrix(0, k + 1, k))
    coefx <- fast.solve(L)%*%m2[, 1]
    coefy <- fast.solve(L)%*%m2[, 2]
    if(k == 3){coefz <- fast.solve(L)%*%m2[, 3]}
    fx <- function(matr, M, coef, step){
      Xn <- numeric(q)
      for (i in 1:q){
        Z <- apply((matr-matrix(M[i,], p, k, byrow = TRUE))^2, 1, sum)
        ifelse(k == 2, Z1<-Z*log(Z), Z1<-sqrt(Z)); Z1[which(is.na(Z1))] <- 0
        ifelse(k == 2, Xn[i] <- coef[p+1] + coef[p+2]*M[i,1] + coef[p+3]*M[i,2] + sum(coef[1:p]*Z1),
               Xn[i] <- coef[p+1] + coef[p+2]*M[i,1] + coef[p+3]*M[i,2] + coef[p+4]*M[i,3] + sum(coef[1:p]*Z1))
        if(PB == TRUE){setTxtProgressBar(pb, step + i)}
      }
      return(Xn)
    }
    matg <- matrix(NA, q, k)
    if(PB==TRUE){pb <- txtProgressBar(min = 0, max = q*k, style = 3) }
    matg[,1] <- fx(matr, M, coefx, step = 1)
    matg[,2] <- fx(matr, M, coefy, step=q)
    if(k==3){matg[,3] <- fx(matr, M, coefz, step=q*2)
    }
    if(PB==TRUE) close(pb)
    return(matg)
  }

  fast.ginv <- function(X, tol = sqrt(.Machine$double.eps)){
    X <- as.matrix(X)
    k <- ncol(X)
    Xsvd <- La.svd(X, k, k)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    rtu <-((1 / Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
    v <-t(Xsvd$vt)[, Positive, drop = FALSE]
    v%*%rtu
  }

  # format
  if (is.list(surface_mesh)) {
    bone <- t(surface_mesh$vb)[, c(1:3)]
  } else {
    bone <- surface_mesh
  }

  # closet vertex to landmark
  lmk.add <- NULL
  for(i in 1:nrow(landmarks)){
    lmk.add <- rbind(lmk.add,
                     which.min(sqrt((landmarks[i, 1] - bone[, 1]) ^ 2 +
                                    (landmarks[i, 2] - bone[, 2]) ^ 2 +
                                    (landmarks[i, 3] - bone[, 3]) ^ 2))[1])
  }

  nlandmarks <- nrow(landmarks)

  # center bone
  bone_centered <- center(bone)
  bone_trans <- colMeans(bone)

  # center template
  template <- center(template) * (csize(bone_centered[lmk.add, ]) / csize(template[(1:nlandmarks), ]))
  template <- template %*% rotate.mat(bone_centered[lmk.add, ], template[(1:nlandmarks), ])

  # sliding points
  template.tps <- tps2d3d(template[-(1:nlandmarks), ], template[(1:nlandmarks), ], bone_centered[lmk.add, ])
  spec.surfs <- bone_centered[-lmk.add, ]
  nei <- numeric(dim(template.tps)[1])
  sliders <- matrix(NA, nrow = dim(template.tps)[1], ncol = 3)
  for (i in 1:dim(template.tps)[1])     {
    nei[i] <- which.min(sqrt((template.tps[i, 1] - spec.surfs[, 1]) ^ 2 +
                             (template.tps[i, 2] - spec.surfs[ ,2]) ^ 2 +
                             (template.tps[i, 3] - spec.surfs[, 3]) ^ 2))[1]
    sliders[i,] <- spec.surfs[nei[i], ]
    spec.surfs <- spec.surfs[-nei[i], ]
  }

  # make output matrix
  selected.out <- rbind(bone_centered[lmk.add, ], sliders)

  # translate back
  selected.out[, 1] <- selected.out[, 1] - (bone_trans[1] * - 1)
  selected.out[, 2] <- selected.out[, 2] - (bone_trans[2] * - 1)
  selected.out[, 3] <- selected.out[, 3] - (bone_trans[3] * - 1)

  # return
  return(selected.out)
}


#' Finds material properties of bone at surface point
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh List. Mesh data imported via ply_import function
#' @param mapped_coords Data frame. 3D coords of remapped surface points
#' @param normal_dist Numeric. Distance surface normal should penetrate surface
#' @param nifti Nifti
#' @param betaCT Numeric. Calibration value for CT to density calculation
#' @param sigmaCT Numeric. Calibration value for CT to density calculation
#' @return Vector. Vector with value for each point on surface
#' @importFrom oro.nifti img_data
#' @importFrom RNifti niftiHeader
#' @export
surface_normal_intersect <- function(surface_mesh, mapped_coords, normal_dist = 3.0, nifti,
                                     betaCT = 1.0, sigmaCT = 1.0) {
  # format surface data
  surface_coords <- t(surface_mesh$vb)[, c(1:3)]
  surface_normals <- t(surface_mesh$normals)[, c(1:3)]

  # format new point data
  vertex_coords <- data.matrix(mapped_coords)
  dims <- dim(vertex_coords)
  vertex_coords <- as.numeric(vertex_coords)
  dim(vertex_coords) <- dims

  # format image data, with voxel coordinates
  #orientation(nifti) <- "RAS"
  img_data <- img_data(nifti)
  dims <- dim(img_data)
  x_seq <- seq(niftiHeader(nifti)$qoffset_x * -1, by = niftiHeader(nifti)$srow_x[1] * -1,
               length.out = dims[1])
  y_seq <- seq(niftiHeader(nifti)$qoffset_y * -1, by = niftiHeader(nifti)$srow_y[2] * -1,
               length.out = dims[2])
  z_seq <- seq(niftiHeader(nifti)$qoffset_z, by = niftiHeader(nifti)$srow_z[3],
               length.out = dims[3])

  # check bone is within scan volume
  bone_x_min <- min(vertex_coords[, 1])
  bone_x_max <- max(vertex_coords[, 1])
  bone_y_min <- min(vertex_coords[, 2])
  bone_y_max <- max(vertex_coords[, 2])
  bone_z_min <- min(vertex_coords[, 3])
  bone_z_max <- max(vertex_coords[, 3])
  vol_x_min <- min(x_seq)
  vol_x_max <- max(x_seq)
  vol_y_min <- min(y_seq)
  vol_y_max <- max(y_seq)
  vol_z_min <- min(z_seq)
  vol_z_max <- max(z_seq)
  x1_good <- bone_x_min > vol_x_min
  x2_good <- bone_x_max < vol_x_max
  y1_good <- bone_y_min > vol_y_min
  y2_good <- bone_y_max < vol_y_max
  z1_good <- bone_z_min > vol_z_min
  z2_good <- bone_z_max < vol_z_max
  vals <- c(x1_good, x2_good, y1_good, y2_good, z1_good, z2_good)
  print(c("x: ", bone_x_min, bone_x_max, vol_x_min, vol_x_max))
  print(c("y: ", bone_y_min, bone_y_max, vol_y_min, vol_y_max))
  print(c("z: ", bone_z_min, bone_z_max, vol_z_min, vol_z_max))
  if (all(vals) != TRUE) {stop("bone not within scan volume")}

  # Find voxels intercepted by line
  mat_peak <- rep(NA, times = nrow(vertex_coords))
  for (i in 1:nrow(vertex_coords)) {
    # find nearest point
    yy <- t(as.matrix(vertex_coords[i, ], 1, 3))
    y <- cdist(yy, surface_coords)
    matched_point <- which.min(y)

    # points to test
    start_point <- surface_coords[matched_point, ]
    end_point <- surface_coords[matched_point, ] + (surface_normals[matched_point, ] * -1 * normal_dist)
    px <- seq(from = start_point[1], to = end_point[1], length.out = 10)
    py <- seq(from = start_point[2], to = end_point[2], length.out = 10)
    pz <- seq(from = start_point[3], to = end_point[3], length.out = 10)

    # for each point along line find voxel with max value
    max_line <- rep(NA, times = 10)
    for (j in 1:10) {
      voxel <- c(which.min(abs(px[j] - x_seq)),
                 which.min(abs(py[j] - y_seq)),
                 which.min(abs(pz[j] - z_seq)))

      max_line[j] <- img_data[voxel[1], voxel[2], voxel[3]]
    }

    # add HU column
    mat_peak[i] <- max(max_line)
  }

  # convert to density
  mat_peak <- ((mat_peak - betaCT) / sigmaCT)

  # return
  return(mat_peak)
}


#' Finds material properties of bone at any point
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param vertex_coords Matrix
#' @param nifti nifti object
#' @param betaCT Calibration value for CT to density calculation
#' @param sigmaCT Calibration value for CT to density calculation
#' @return Vector. Vector with value for each point on surface
#' @export
voxel_point_intersect <- function(vertex_coords, nifti, betaCT, sigmaCT) {
  vertex_coords <- data.matrix(vertex_coords)
  dims <- dim(vertex_coords)
  vertex_coords <- as.numeric(vertex_coords)
  dim(vertex_coords) <- dims

  # format image data, with voxel coordinates
  img_data <- img_data(nifti)
  dims <- dim(img_data)
  x_seq <- seq(niftiHeader(nifti)$qoffset_x * -1, by = niftiHeader(nifti)$srow_x[1] * -1,
               length.out = dims[1])
  y_seq <- seq(niftiHeader(nifti)$qoffset_y * -1, by = niftiHeader(nifti)$srow_y[2] * -1,
               length.out = dims[2])
  z_seq <- seq(niftiHeader(nifti)$qoffset_z, by = niftiHeader(nifti)$srow_z[3],
               length.out = dims[3])

  ## check vertices are in scan volume
  bone_scan_check(vertex_coords, nifti)

  ## Find voxels intercepted by line
  mat_peak <- rep(NA, times = nrow(vertex_coords))
  for (i in 1:nrow(vertex_coords)) {
    # start and end points of line
    point <- vertex_coords[i, ]

    # find voxel
    voxel <- c(which.min(abs(point[1] - x_seq)),
               which.min(abs(point[2] - y_seq)),
               which.min(abs(point[3] - z_seq)))

    # add HU column
    mat_peak[i] <- img_data[voxel[1], voxel[2], voxel[3]]
  }

  # convert to density
  mat_peak <- ((mat_peak - betaCT) / sigmaCT)

  # return
  return(mat_peak)
}


#' Checks bone is in scan volume
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param vertex_coords Matrix. 3D bone vertex coordinates
#' @param nifti nifti image
#' @return Error message if not in scan volume
#' @export
bone_scan_check <- function(vertex_coords, nifti) {
  # format image data, with voxel coordinates
  img_data <- img_data(nifti)
  dims <- dim(img_data)
  x_seq <- seq(niftiHeader(nifti)$qoffset_x * -1, by = niftiHeader(nifti)$srow_x[1] * -1,
               length.out = dims[1])
  y_seq <- seq(niftiHeader(nifti)$qoffset_y * -1, by = niftiHeader(nifti)$srow_y[2] * -1,
               length.out = dims[2])
  z_seq <- seq(niftiHeader(nifti)$qoffset_z, by = niftiHeader(nifti)$srow_z[3],
               length.out = dims[3])

  # check bone is within scan volume
  bone_x_min <- min(vertex_coords[, 1])
  bone_x_max <- max(vertex_coords[, 1])
  bone_y_min <- min(vertex_coords[, 2])
  bone_y_max <- max(vertex_coords[, 2])
  bone_z_min <- min(vertex_coords[, 3])
  bone_z_max <- max(vertex_coords[, 3])
  vol_x_min <- min(x_seq)
  vol_x_max <- max(x_seq)
  vol_y_min <- min(y_seq)
  vol_y_max <- max(y_seq)
  vol_z_min <- min(z_seq)
  vol_z_max <- max(z_seq)
  x1_good <- bone_x_min > vol_x_min
  x2_good <- bone_x_max < vol_x_max
  y1_good <- bone_y_min > vol_y_min
  y2_good <- bone_y_max < vol_y_max
  z1_good <- bone_z_min > vol_z_min
  z2_good <- bone_z_max < vol_z_max
  vals <- c(x1_good, x2_good, y1_good, y2_good, z1_good, z2_good)
  print(c("x: ", bone_x_min, bone_x_max, vol_x_min, vol_x_max))
  print(c("y: ", bone_y_min, bone_y_max, vol_y_min, vol_y_max))
  print(c("z: ", bone_z_min, bone_z_max, vol_z_min, vol_z_max))
  if (all(vals) != TRUE) {stop("bone not within scan volume")}
}


#' Finds point closest to vertex for all vertices in a surface mesh
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh mesh object
#' @param template_points matrix
#' @return Vector. Closest point on mesh to each template pint
#' @importFrom rdist cdist
#' @export
mesh_template_match <- function(surface_mesh, template_points) {
  # get vertex coords
  vertex_coords <- t(surface_mesh$vb)[, c(1:3)]

  # calculate distances
  y <- cdist(vertex_coords, template_points)

  # calculate which point is closest
  matched_points <- apply(y, 1, which.min)

  # return
  return(matched_points)
}


#' maps numeric values to a color
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param x Vector.
#' @param maxi Numeric.
#' @param mini Numeric.
#' @return Vector of same length as x
#' @importFrom grDevices colorRamp rgb
#' @export
color_mapping <- function(x, maxi, mini) {
  # scale distance vector between 0 and 1
  if (missing(maxi) == TRUE & missing(mini) == TRUE) {
    x01 <- (x - min(x)) / (max(x) - min(x))
  }
  if (missing(maxi) == FALSE & missing(mini) == FALSE) {
    x01 <- (x - mini) / (maxi - mini)
  }
  x01[is.na(x01)] <- mini

  # color map
  colormap <- colorRamp(c("dark blue", "blue", "light blue",
                          "green", "yellow", "red", "pink"),
                        interpolate = "spline")
  ply_col <- colormap(x01)
  ply_col <- apply(ply_col, 1, function(x) rgb(x[1], x[2], x[3],
                                               maxColorValue = 255))

  # return colour map frame
  return(ply_col)
}


#' plot mesh
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh Mesh object
#' @param density_color Vector
#' @param title String
#' @param userMat Matrix
#' @return plot of mesh with color
#' @importFrom rgl shade3d view3d bgplot3d
#' @importFrom graphics plot.new mtext
#' @importFrom methods hasArg
#' @export
plot_mesh <- function(surface_mesh, density_color, title, userMat) {
  verts <- as.matrix(t(surface_mesh$vb)[,-4])
  surface_mesh$vb <- rbind(t(verts), 1)
  surface_mesh$material <- list(color = density_color)
  shade3d(surface_mesh, meshColor = "vertices", main = "Age", specular = 'black')
  bgplot3d({
    plot.new()
    mtext(side = 1, title, line = 3)
  })
  if (hasArg(userMat)) {
    view3d(userMatrix = userMat)
  }
}


#' Produce stand alone color bar
#' @param colors String
#' @param mini Numeric
#' @param maxi Numeric
#' @param breaks Numeric vector
#' @param orientation "horizontal" or "vertical"
#' @param title String
#' @param text_size Numeric
#' @param plot Logical
#' @importFrom ggplot2 ggplot unit labs guides theme element_text geom_point aes scale_color_gradientn guide_colorbar
#' @importFrom cowplot get_legend
#' @importFrom ggpubr as_ggplot
#' @export
color_bar <- function(colors, mini, maxi, orientation = "vertical", breaks,
                      title = "", text_size = 11, plot = TRUE) {
  # packages needed
  #require(ggplot2)
  #require(cowplot)
  #require(ggpubr)
  x <- y <- NULL

  # make plot
  z2 <- seq(from = mini, to = maxi, length.out = 100)
  df <- data.frame(x = 1:100, y = 1:100, z2 = z2)
  g <- ggplot(df, aes(x, y)) + geom_point(aes(color = z2))
  g <- g + scale_color_gradientn(colors = colors,
                                 breaks = breaks)
  g <- g + labs(color = title)
  if (orientation == "horizontal") {
    g <- g + guides(color = guide_colorbar(title.position = "top"))
    g <- g + theme(legend.key.size = unit(2.0, "cm"),
                   legend.position = "bottom")
  }
  if (orientation == "vertical") {
    g <- g + theme(legend.key.size = unit(2.0, "cm"),
                   legend.text = element_text(size = text_size))
  }

  # extract legend
  legend <- get_legend(g)
  lg <- as_ggplot(legend)

  # plot if requested
  if (plot == TRUE) {
    lg
  }

  # return legend
  return(lg)
}


#' fill bone
#' @param bone_surface Mesh object
#' @param spacing Numeric
#' @return Matrix with internal point coordinates
#' @importFrom ptinpoly pip3d
#' @importFrom stats runif
#' @export
fill_bone_points <- function(bone_surface, spacing) {
  # verts
  verts <- t(bone_surface$vb)[, 1:3]
  faces <- t(bone_surface$it)

  # bone extents
  x_min <- min(verts[, 1])
  x_max <- max(verts[, 1])
  y_min <- min(verts[, 2])
  y_max <- max(verts[, 2])
  z_min <- min(verts[, 3])
  z_max <- max(verts[, 3])

  # make point df
  x <- seq(from = x_min, to = x_max, by = spacing)
  y <- seq(from = y_min, to = y_max, by = spacing)
  z <- seq(from = z_min, to = z_max, by = spacing)
  pt_mat <- as.matrix(expand.grid(x = x, y = y, z = z))

  # find which points are within surface
  in_bone <- which(pip3d(verts, faces, pt_mat) == 1)
  in_coords <- pt_mat[in_bone, ]

  # add a little noise
  dims <- dim(in_coords)
  noise <- runif(dims[1] * dims[2], -0.001, 0.001)
  noise_mat <- matrix(noise, nrow = dims[1], ncol = dims[2])
  in_coords <- in_coords + noise_mat

  # return
  return(in_coords)
}


#' Color mesh
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh Mesh object
#' @param template_pts Matrix
#' @param density_vector Vector
#' @param maxi Numeric
#' @param mini Numeric
#' @param export_path Character
#' @importFrom Rvcg vcgPlyWrite
#' @export
color_mesh <- function(surface_mesh, template_pts, density_vector, maxi = 2000,
                       mini = 0, export_path) {
  # mesh match
  mesh_match <- mesh_template_match(surface_mesh, template_pts)

  # color
  density_vector[density_vector > maxi] <- maxi
  density_vector[density_vector < mini] <- mini
  color_map <- color_mapping(density_vector, maxi, mini)

  # color to mesh
  surface_mesh$material$color <- color_map[mesh_match]

  # return
  if (missing(export_path) == FALSE) {
    vcgPlyWrite(surface_mesh, export_path)
  } else {
    return(surface_mesh)
  }
}


#' local significance
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param vertices Matrix
#' @param sig_vals Numeric vector
#' @param sig_level Numeric. Default 0.05
#' @param dist Numeric. Distance to check for vertices
#' @param n_local Numeric. Number of local significant values needed
#' @return Numeric vector
#' @importFrom rdist cdist
#' @export
rm_local_sig <- function(vertices, sig_vals, sig_level = 0.05, dist, n_local = 1) {
  # identify significant values
  sig_inds <- which(sig_vals < sig_level)

  # check if nearby values are also significant
  sig_vals_updated <- sig_vals
  for (i in seq_along(sig_inds)) {
    # which vertices are within distance
    vert <- vertices[sig_inds[i], ]
    y <- cdist(vert, vertices)
    z <- which(y < dist)

    # are these also significant? If so update vector
    #print(sum(sig_vals[z] < sig_level))
    if (sum(sig_vals[z] < sig_level) < n_local) {sig_vals_updated[i] = sig_level + 0.01}
  }

  # return vector
  return(sig_vals_updated)
}

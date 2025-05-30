# To-do
# add .nrrd compatibility
# add multiple bone scan/mesh/landmark sets
# optimize surface_normal_intersect

#' import landmark coordinates
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param landmark_path String. File path to landmark data, .json or .fcsv format
#' @return dataframe. Columns are landmark name, x, y, and z coordinates
#' @examples
#' landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                              package = "BoneDensityMapping")
#' landmarks <- import_lmks(landmark_path)
#' @importFrom rjson fromJSON
#' @importFrom tools file_ext
#' @export
import_lmks <- function(landmark_path) {
  file_type <- file_ext(landmark_path)

  if (file_type == "json") {
    lmks <- fromJSON(file = landmark_path)

    # extract point lists
    lmks_ <- lmks$markups[[1]]$controlPoints

    # extract names and positions
    lmk_names <- rep(NA, length(lmks_))
    coords <- matrix(NA, nrow = length(lmks_), ncol = 3)
    for (i in seq_along(lmks_)) {
      fid <- lmks_[[i]]
      lmk_names[i] <- fid$label
      coords[i, ] <- fid$position
    }
    df <- data.frame(lmk_id = lmk_names, x = coords[, 1], y = coords[, 2], z = coords[, 3])

  } else if (file_type == "fcsv") {
    coords <- read.csv(landmark_path, skip = 3, header = FALSE)[, 2:4]
    df <- data.frame(lmk_id = seq_len(nrow(coords)), x = coords[[1]], y = coords[[2]], z = coords[[3]])

  } else {
    stop("Unsupported file type: must be .json or .fcsv")
  }

  # return
  return(df)
}


#' import CT scan
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param scan_filepath String. File path to CT scan data. Should be .nii or .nrrd
#' @return scan object
#' @examples
#' scan_filepath <- system.file("extdata", "test_CT_hip.nii",
#'                          package = "BoneDensityMapping")
#' import_scan(scan_filepath)
#' @importFrom oro.nifti readNIfTI
#' @export
import_scan <- function(scan_filepath) {
  file_type <- file_ext(scan_filepath)

  if (file_type == "nii") {
    nifti <- readNIfTI(scan_filepath)
  }
  else {
    stop("Unsupported file type: must be .nii")
  }
  return(nifti)
}


#' import surface mesh
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh_filepath String. File path to CT scan data. Should be .stl or .ply
#' @return mesh object
#' @examples
#' surface_mesh_filepath <- system.file("extdata", "test_CT_femur.stl",
#'                          package = "BoneDensityMapping")
#' surface_mesh <- import_mesh(surface_mesh_filepath)
#' @importFrom Rvcg vcgImport
#' @export
import_mesh <- function(surface_mesh_filepath) {
  # import scan
  surface_mesh <- vcgImport(surface_mesh_filepath)

  # return
  return(surface_mesh)
}


#' check landmarks are close to the mesh
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh mesh object
#' @param landmarks dataframe. Columns are landmark name, x, y, and z coordinates
#' @param threshold Numeric. Distance landmark can be from surface without
#' warning being thrown
#' @return String. Returns a message warning that landmarks are not on bone
#' surface
#' @examples
#' landmark_path <- system.file("extdata", "test_femur.fcsv",
#'                              package = "BoneDensityMapping")
#' landmarks <- import_lmks(landmark_path)
#' surface_mesh_filepath <- system.file("extdata", "test_CT_femur.stl",
#'                                  package = "BoneDensityMapping")
#' surface_mesh <- import_mesh(surface_mesh_filepath)
#' landmark_check(surface_mesh, landmarks, threshold = 1.0)
#' @importFrom rdist cdist
#' @importFrom utils read.csv
#' @export
landmark_check <- function(surface_mesh, landmarks, threshold = 1.0) {
  vertices <- t(surface_mesh$vb)[, c(1:3)]
  coords <- landmarks[, c("x", "y", "z")]

  dists <- c()
  for (i in 1:nrow(coords)) {
    pt <- matrix(as.numeric(unlist(coords[i, ])), nrow = 1)
    x <- cdist(vertices, pt)
    dists <- c(dists, min(x))
  }

  # return message if landmarks not on bone surface
  if (any(dists > threshold)) {
    bad_ids <- landmarks$lmk_id[dists > threshold]
    message("landmarks not on bone surface: ", paste(bad_ids, collapse = ", "))
  }
  else message("landmarks are on bone surface")
}


#' Check if mesh or vertex coordinates are fully contained within scan volume
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh A 3D mesh object (class \code{mesh3d}) or a numeric
#'    matrix/data frame of vertex coordinates (with columns representing X, Y, and Z).
#' @param nifti A NIfTI image object representing a CT scan.
#' @return If all vertices are inside the scan volume, the function prints a summary
#'   of the mesh and scan bounds. If any vertices lie outside the scan volume, it
#'   raises an error.
#' @examples
#' scan_filepath <- system.file("extdata", "test_CT_hip.nii",
#'                           package = "BoneDensityMapping")
#' nifti <- import_scan(scan_filepath)
#' surface_mesh_filepath <- system.file("extdata", "test_CT_femur.stl",
#'                                  package = "BoneDensityMapping")
#' surface_mesh <- import_mesh(surface_mesh_filepath)
#' bone_scan_check(surface_mesh, nifti)
#' @export
bone_scan_check <- function(surface_mesh, nifti) {

  if (inherits(surface_mesh, "mesh3d")) {
    vertices <- t(surface_mesh$vb)[, 1:3]
  } else if (is.matrix(surface_mesh) || is.data.frame(surface_mesh)) {
    vertices <- as.matrix(surface_mesh)
  } else {
    stop("surface_mesh must be a mesh3d object or a matrix of vertex coordinates.")
  }

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
  mesh_x_min <- min(vertices[, 1])
  mesh_x_max <- max(vertices[, 1])
  mesh_y_min <- min(vertices[, 2])
  mesh_y_max <- max(vertices[, 2])
  mesh_z_min <- min(vertices[, 3])
  mesh_z_max <- max(vertices[, 3])
  vol_x_min <- min(x_seq)
  vol_x_max <- max(x_seq)
  vol_y_min <- min(y_seq)
  vol_y_max <- max(y_seq)
  vol_z_min <- min(z_seq)
  vol_z_max <- max(z_seq)
  x1_good <- mesh_x_min > vol_x_min
  x2_good <- mesh_x_max < vol_x_max
  y1_good <- mesh_y_min > vol_y_min
  y2_good <- mesh_y_max < vol_y_max
  z1_good <- mesh_z_min > vol_z_min
  z2_good <- mesh_z_max < vol_z_max
  vals <- c(x1_good, x2_good, y1_good, y2_good, z1_good, z2_good)
  df <- data.frame(
    Axis = c("X", "Y", "Z"),
    Mesh_Min = c(mesh_x_min, mesh_y_min, mesh_z_min),
    Mesh_Max = c(mesh_x_max, mesh_y_max, mesh_z_max),
    Scan_Min = c(vol_x_min, vol_y_min, vol_z_min),
    Scan_Max = c(vol_x_max, vol_y_max, vol_z_max)
  )
  print(df, digits = 3)
  if (!all(vals)) {stop("Mesh not within scan volume")}
}


#' fill bone
#' @param surface_mesh Mesh object
#' @param spacing Numeric
#' @return Matrix with internal point coordinates
#' @examples
#' surface_mesh_filepath <- system.file("extdata", "test_CT_femur.stl",
#'                                  package = "BoneDensityMapping")
#' surface_mesh <- import_mesh(surface_mesh_filepath)
#' internal_fill <- fill_bone_points(surface_mesh, 10)
#' @importFrom ptinpoly pip3d
#' @importFrom stats runif
#' @export
fill_bone_points <- function(surface_mesh, spacing) {
  # verts
  vertices <- t(surface_mesh$vb)[, 1:3]
  faces <- t(surface_mesh$it)

  # bone extents
  x_min <- min(vertices[, 1])
  x_max <- max(vertices[, 1])
  y_min <- min(vertices[, 2])
  y_max <- max(vertices[, 2])
  z_min <- min(vertices[, 3])
  z_max <- max(vertices[, 3])

  # make point df
  x <- seq(from = x_min, to = x_max, by = spacing)
  y <- seq(from = y_min, to = y_max, by = spacing)
  z <- seq(from = z_min, to = z_max, by = spacing)
  pt_mat <- as.matrix(expand.grid(x = x, y = y, z = z))

  # find which points are within surface
  in_bone <- which(pip3d(vertices, faces, pt_mat) == 1)
  in_coords <- pt_mat[in_bone, ]

  # add a little noise
  dims <- dim(in_coords)
  noise <- runif(dims[1] * dims[2], -0.001, 0.001)
  noise_mat <- matrix(noise, nrow = dims[1], ncol = dims[2])
  in_coords <- in_coords + noise_mat

  # return
  return(in_coords)
}

#' Redefine surface points. Adds additional surface points (“sliders”) that are spatially distributed across the mesh surface.
#' @author Scott Telfer \email{scott.telfer@gmail.com} Adapted from geomorph
#' @param surface_mesh Mesh object
#' @param landmarks Data frame with landmark coordinates (columns: ID, x, y, z)
#' @param no_surface_sliders Numeric, number of additional surface points to generate
#' @return Data frame. 3D coordinates for the combined set of original landmarks and the new surface points
#' @examples
#' surface_mesh_path <- system.file("extdata", "test_CT_femur.stl",
#'                                  package = "BoneDensityMapping")
#' surface_mesh <- import_mesh(surface_mesh_path)
#' landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                              package = "BoneDensityMapping")
#' landmarks <- import_lmks(landmark_path)
#' mapped_coords <- surface_points_template(surface_mesh, landmarks, 1000)
#' @importFrom stats kmeans
#' @export
surface_points_template <- function(surface_mesh, landmarks, no_surface_sliders) {
  # extract vertex coordinates from mesh
  if (is.list(surface_mesh)) {
    vertices <- t(surface_mesh$vb)[, c(1:3)]
  } else {
    vertices <- surface_mesh
  }
  colnames(vertices) <- c("xpts", "ypts", "zpts")

  # Isolate just the x, y, z coordinates from the landmark df
  landmark_coords <- landmarks[, 2:4]
  landmark_coords <- apply(landmark_coords, 2, as.numeric)  # ensure numeric

  # Identify which mesh vertices are closest to each landmark, add them to lmk_add, remove these mesh vertices
  lmk.add <- NULL
  for(i in 1:nrow(landmark_coords)){
    lmk.add <- rbind(lmk.add, which.min(sqrt((landmark_coords[i, 1] - vertices[, 1]) ^ 2 +
                                               (landmark_coords[i, 2] - vertices[, 2]) ^ 2 +
                                               (landmark_coords[i, 3] - vertices[, 3]) ^ 2))[1])
  }
  vertices <- vertices[-lmk.add, ]

  # Use k-means clustering to find 'no_surface_sliders' new surface points
  colnames(landmark_coords) <- c("xpts", "ypts", "zpts")
  new_surface_points <- rbind(landmark_coords,
                              kmeans(x = vertices, centers = no_surface_sliders,
                                     iter.max = 25)$centers)

  # return
  return(as.data.frame(new_surface_points))
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
      if (inherits(res, "try-error")) res <- fast.ginv(x)
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
#' @param surface_mesh Mesh object
#' @param mapped_coords Data frame. 3D coords of remapped surface points
#' @param normal_dist Numeric. Distance surface normal should penetrate surface
#' @param nifti Nifti CT scan image
#' @param betaCT Numeric. Calibration value for CT to density calculation
#' @param sigmaCT Numeric. Calibration value for CT to density calculation
#' @param rev_x Logical.
#' @param rev_y Logical.
#' @param rev_z Logical.
#' @return Vector. Vector with value for each point on surface
#' @examples
#' surface_mesh_path <- system.file("extdata", "test_CT_femur.stl",
#'                                  package = "BoneDensityMapping")
#' surface_mesh <- import_mesh(surface_mesh_path)
#' landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                              package = "BoneDensityMapping")
#' landmarks <- import_lmks(landmark_path)
#' nifti_path <- system.file("extdata", "test_CT_hip.nii",
#'                           package = "BoneDensityMapping")
#' nifti <- import_scan(nifti_path)
#' mapped_coords <- surface_points_template(surface_mesh, landmarks, no_surface_sliders = 100)
#' mat_peak <- surface_normal_intersect(surface_mesh, mapped_coords,
#'                                      normal_dist = 3.0, nifti, betaCT = 1.0,
#'                                      sigmaCT = 1.0)
#' @importFrom oro.nifti img_data
#' @importFrom RNifti niftiHeader
#' @export
surface_normal_intersect <- function(surface_mesh, mapped_coords, normal_dist = 3.0,
                                     nifti, betaCT = 1.0, sigmaCT = 1.0,
                                     rev_x = FALSE, rev_y = FALSE, rev_z = FALSE) {
  surface_coords <- t(surface_mesh$vb)[, c(1:3)]
  surface_normals <- t(surface_mesh$normals)[, c(1:3)]

  # Convert mapped_coords to numeric matrix form for calculations
  vertices <- data.matrix(mapped_coords)
  dims <- dim(vertices)
  vertices <- as.numeric(vertices)
  dim(vertices) <- dims

  # format image data, with voxel coordinates
  img_data <- img_data(nifti)
  dims <- dim(img_data)
  x_by <- niftiHeader(nifti)$srow_x[1]
  y_by <- niftiHeader(nifti)$srow_y[2]
  z_by <- niftiHeader(nifti)$srow_z[3]
  if (rev_x == TRUE) {
    x_seq <- rev(seq(niftiHeader(nifti)$qoffset_x * -1, by = x_by * -1, length.out = dims[1]))
  } else {
    x_seq <- seq(niftiHeader(nifti)$qoffset_x * -1, by = x_by * -1, length.out = dims[1])
  }
  if (rev_y == TRUE) {
    y_seq <- rev(seq(niftiHeader(nifti)$qoffset_y * -1, by = y_by * -1, length.out = dims[2]))
  } else {
    y_seq <- seq(niftiHeader(nifti)$qoffset_y * -1, by = y_by * -1, length.out = dims[2])
  }
  if (rev_z == TRUE) {
    z_seq <- rev(seq(niftiHeader(nifti)$qoffset_z, by = z_by, length.out = dims[3]))
  } else {
    z_seq <- seq(niftiHeader(nifti)$qoffset_z, by = z_by, length.out = dims[3])
  }

  # check bone is within scan volume
  mesh_x_min <- min(vertices[, 1])
  mesh_x_max <- max(vertices[, 1])
  mesh_y_min <- min(vertices[, 2])
  mesh_y_max <- max(vertices[, 2])
  mesh_z_min <- min(vertices[, 3])
  mesh_z_max <- max(vertices[, 3])
  vol_x_min <- min(x_seq)
  vol_x_max <- max(x_seq)
  vol_y_min <- min(y_seq)
  vol_y_max <- max(y_seq)
  vol_z_min <- min(z_seq)
  vol_z_max <- max(z_seq)
  x1_good <- mesh_x_min > vol_x_min
  x2_good <- mesh_x_max < vol_x_max
  y1_good <- mesh_y_min > vol_y_min
  y2_good <- mesh_y_max < vol_y_max
  z1_good <- mesh_z_min > vol_z_min
  z2_good <- mesh_z_max < vol_z_max
  vals <- c(x1_good, x2_good, y1_good, y2_good, z1_good, z2_good)
  df <- data.frame(
    Axis = c("X", "Y", "Z"),
    Mesh_Min = c(mesh_x_min, mesh_y_min, mesh_z_min),
    Mesh_Max = c(mesh_x_max, mesh_y_max, mesh_z_max),
    Scan_Min = c(vol_x_min, vol_y_min, vol_z_min),
    Scan_Max = c(vol_x_max, vol_y_max, vol_z_max)
  )
  print(df, digits = 3)
  if (!all(vals)) {stop("Mesh not within scan volume")}

  # Find voxels intercepted by line
  mat_peak <- rep(NA, times = nrow(vertices))

  # Loop through each surface point to sample CT data along the surface normal
  for (i in 1:nrow(vertices)) {
    # find nearest point
    yy <- t(as.matrix(vertices[i, ], 1, 3))
    y <- cdist(yy, surface_coords)
    matched_point <- which.min(y)

    # Define start and end points along the surface normal line (penetrate normal_dist mm inside).
    # Interpolate 10 points along this line inside the bone surface
    start_point <- surface_coords[matched_point, ]
    end_point <- surface_coords[matched_point, ] + (surface_normals[matched_point, ] * -1 * normal_dist)
    px <- seq(from = start_point[1], to = end_point[1], length.out = 10)
    py <- seq(from = start_point[2], to = end_point[2], length.out = 10)
    pz <- seq(from = start_point[3], to = end_point[3], length.out = 10)

    # For each point along the line, find corresponding voxel and record CT intensity
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

  # Convert Hounsfield units (HU) to density using calibration values
  mat_peak <- ((mat_peak - betaCT) / sigmaCT)

  # Return the vector of density values for each mapped coordinate
  return(mat_peak)
}


#' Finds material properties of bone at any point
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param vertex_coords Matrix
#' @param nifti nifti object
#' @param betaCT Calibration value for CT to density calculation
#' @param sigmaCT Calibration value for CT to density calculation
#' @param check_in_vol Logical Include check that model is in scans volume
#' and print dimensions
#' @examples
#' surface_mesh_path <- system.file("extdata", "test_CT_femur.stl",
#'                                  package = "BoneDensityMapping")
#' surface_mesh <- import_mesh(surface_mesh_path)
#' landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                              package = "BoneDensityMapping")
#' landmarks <- import_lmks(landmark_path)
#' nifti_path <- system.file("extdata", "test_CT_hip.nii",
#'                           package = "BoneDensityMapping")
#' nifti <- import_scan(nifti_path)
#' mapped_coords <- surface_points_template(surface_mesh, landmarks, no_surface_sliders = 100)
#' mat_peak <- voxel_point_intersect(mapped_coords, nifti,
#'    betaCT = 1.0, sigmaCT = 1.0, check_in_vol = TRUE)
#' @return Vector. Vector with value for each point on surface
#' @export
voxel_point_intersect <- function(vertex_coords, nifti, betaCT = 1.0, sigmaCT = 1.0,
                                  check_in_vol = FALSE) {
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
  if (check_in_vol == TRUE) {
    bone_scan_check(vertex_coords, nifti)
  }

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


#' Finds the closest template point for each vertex in a 3D surface mesh
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh mesh object
#' @param template_points matrix
#' @return Vector. Closest point on mesh to each template point
#' @examples
#' surface_mesh_path <- system.file("extdata", "test_CT_femur.stl",
#'                                  package = "BoneDensityMapping")
#' surface_mesh <- import_mesh(surface_mesh_path)
#' landmark_path <- system.file("extdata", "test_femur.mrk.json", package = "BoneDensityMapping")
#' landmarks <- import_lmks(landmark_path)
#' mapped_coords <- surface_points_template(surface_mesh, landmarks, 100)
#' matched_points <- mesh_template_match(surface_mesh, mapped_coords)
#' @importFrom FNN get.knnx
#' @export
mesh_template_match <- function(surface_mesh, template_points) {
  # Get vertex coordinates from the mesh
  vertex_coords <- t(surface_mesh$vb)[, c(1:3)]

  # Ensure template_points is a matrix (not a data frame)
  template_points <- as.matrix(template_points)

  # Use get.knnx to find nearest neighbors
  # This finds the nearest (k = 1) template point to each vertex
  nn <- get.knnx(data = template_points, query = vertex_coords, k = 1)

  # nn$nn.index contains the indices of the nearest neighbors
  return(as.vector(nn$nn.index))
}


#' maps numeric values to a color
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param x Vector.
#' @param maxi Numeric.
#' @param mini Numeric.
#' @param color_sel Vector.
#' @examples
#' surface_mesh_filepath <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")
#' surface_mesh <- import_mesh(surface_mesh_filepath)
#' vertices <- t(surface_mesh$vb)[, c(1:3)]
#' nifti_path <- system.file("extdata", "test_CT_hip.nii", package = "BoneDensityMapping")
#' nifti <- import_scan(nifti_path)
#' mat_peak <- voxel_point_intersect(vertices, nifti)
#' colors <- color_mapping(mat_peak)
#' @return Vector of same length as x
#' @importFrom grDevices colorRamp rgb
#' @export
color_mapping <- function(x, maxi, mini, color_sel) {
  # Scale input vector x between 0 and 1
  if (missing(maxi) & missing(mini)) {
    x01 <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  } else if (!missing(maxi) & !missing(mini)) {
    x01 <- (x - mini) / (maxi - mini)
  } else {
    stop("Either provide both 'maxi' and 'mini' or neither.")
  }

  # Replace NA values with 0
  x01[is.na(x01)] <- 0

  # Clamp values between 0 and 1 to avoid colorRamp errors
  x01[x01 < 0] <- 0
  x01[x01 > 1] <- 1

  # Generate color map
  if (missing(color_sel)) {
    colormap <- colorRamp(c("dark blue", "blue", "light blue",
                            "green", "yellow", "red", "pink"),
                          interpolate = "spline")
  } else {
    colormap <- colorRamp(color_sel)
  }

  # Apply color map to scaled values
  ply_col <- colormap(x01)
  ply_col <- apply(ply_col, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))

  return(ply_col)
}


#' Takes matrix of points (template_points) of length n with
#' density vector of length n and maps it to a surface mesh of length m.
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh Mesh object
#' @param template_pts Matrix
#' @param density_vector Vector
#' @param maxi Numeric
#' @param mini Numeric
#' @param export_path Character
#' @return Mesh with added color dimension
#' @examples
#' surface_mesh_path <- system.file("extdata", "test_CT_femur.stl",
#'                                  package = "BoneDensityMapping")
#' surface_mesh <- import_mesh(surface_mesh_path)
#' landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                              package = "BoneDensityMapping")
#' landmarks <- import_lmks(landmark_path)
#' nifti_path <- system.file("extdata", "test_CT_hip.nii",
#'                           package = "BoneDensityMapping")
#' nifti <- import_scan(nifti_path)
#' mapped_coords <- surface_points_template(surface_mesh, landmarks, no_surface_sliders = 100)
#' mat_peak <- voxel_point_intersect(mapped_coords, nifti, betaCT = 1.0, sigmaCT = 1.0)
#' colored_mesh <- color_mesh(surface_mesh, mapped_coords, mat_peak, maxi=2000, mini=0)
#' @importFrom Rvcg vcgPlyWrite
#' @export
color_mesh <- function(surface_mesh, template_pts, density_vector, maxi = 2000,
                       mini = 0, export_path) {
  # this function
  # then use plot<-color_mesh(surface_mesh, mapped_coords, mat_peak, maxi=2000, mini=0)
  # shade3d(plot)
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


#' plot mesh
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param surface_mesh Mesh object
#' @param density_color Vector. Colors mapped from density values.
#' @param title String. Plot title.
#' @param userMat Optional matrix. Controls graph orientation.
#' @return plot of mesh with color
#' @examples
#' surface_mesh_path <- system.file("extdata", "test_CT_femur.stl",
#'                                  package = "BoneDensityMapping")
#' surface_mesh <- import_mesh(surface_mesh_path)
#' landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                              package = "BoneDensityMapping")
#' landmarks <- import_lmks(landmark_path)
#' nifti_path <- system.file("extdata", "test_CT_hip.nii",
#'                           package = "BoneDensityMapping")
#' nifti <- import_scan(nifti_path)
#' mapped_coords <- surface_points_template(surface_mesh, landmarks, no_surface_sliders = 100)
#' mat_peak <- voxel_point_intersect(mapped_coords, nifti, betaCT = 1.0, sigmaCT = 1.0)
#' colored_mesh <- color_mesh(surface_mesh, mapped_coords, mat_peak, maxi=2000, mini=0)
#' plot <- plot_mesh(colored_mesh, "Femur Bone Density")
#' @importFrom rgl shade3d view3d bgplot3d
#' @importFrom graphics plot.new mtext
#' @importFrom methods hasArg
#' @export
plot_mesh <- function(surface_mesh, title, density_color = NULL, userMat = NULL) {
  #ensure formatting
  vertices <- as.matrix(t(surface_mesh$vb)[,-4])
  surface_mesh$vb <- rbind(t(vertices), 1)

  if (!missing(density_color)) {
    surface_mesh$material <- list(color = density_color)
    shade3d(surface_mesh, meshColor = "vertices", main = "Age", specular = 'black')
  } else {
    shade3d(surface_mesh)
  }

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
#' @examples
#' colors <- c("darkblue", "blue", "lightblue", "green", "yellow", "red", "pink")
#' color_bar(colors, 0, 1, breaks = c(0, 0.25, 0.5, 0.75, 1))
#' @importFrom ggplot2 ggplot unit labs guides theme element_text geom_point aes scale_color_gradientn guide_colorbar
#' @importFrom cowplot get_legend
#' @importFrom ggpubr as_ggplot
#' @export
color_bar <- function(colors, mini, maxi, orientation = "vertical", breaks,
                      title = "", text_size = 11, plot = TRUE) {
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

#' local significance
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param vertices Matrix
#' @param sig_vals Numeric vector
#' @param changes Numeric vector
#' @param sig_level Numeric. Default 0.05
#' @param dist Numeric. Distance to check for vertices
#' @return Numeric vector
#' @importFrom rdist cdist
rm_local_sig <- function(vertices, sig_vals, changes, sig_level = 0.05, dist) {
  # identify significant values
  sig_inds <- which(sig_vals < sig_level)
  sig_changes <- changes[sig_inds]
  sig_inds_up <- sig_inds[which(sig_changes > 0)]
  sig_inds_down <- sig_inds[which(sig_changes < 0)]

  # check if nearby values are also significant
  sig_vals_updated <- sig_vals
  sig_verts_up <- vertices[sig_inds_up, ]
  sig_verts_down <- vertices[sig_inds_down, ]
  for (i in 1:length(sig_inds_up)) {
    vert <- vertices[sig_inds_up[i], ]
    y <- cdist(vert, sig_verts_up)
    if (length(which(y < dist)) < 2) {sig_vals_updated[sig_inds_up[i]] = 0.1}
  }
  for (i in 1:length(sig_inds_down)) {
    vert <- vertices[sig_inds_down[i], ]
    y <- cdist(vert, sig_verts_down)
    if (length(which(y < dist)) < 2) {sig_vals_updated[sig_inds_down[i]] = 0.1}
  }

  # return vector
  return(sig_vals_updated)
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

#' reorientate_landmarks
#' @author Scott Telfer \email{scott.telfer@gmail.com}
#' @param landmark_path String
#' @param x Integer Value to apply to convert mesh i.e. -1 will mirror x coords
#' @param y Integer Value to apply to convert mesh i.e. -1 will mirror y coords
#' @param z Integer Value to apply to convert mesh i.e. -1 will mirror z coords
#' @return Overwritten landmark file
#' @examples
#' landmark_path <- system.file("extdata", "test_femur.mrk.json",
#'                              package = "BoneDensityMapping")
#' reoriented_landmarks <- reorientate_landmarks(landmark_path)
#' @importFrom utils read.table write.table
#' @importFrom jsonlite read_json write_json
#' @export
reorientate_landmarks <- function(landmark_path, x = 1, y = 1, z = 1) {
  # Check file extension
  file_ext <- tools::file_ext(landmark_path)

  if (file_ext == "fcsv") {
    # Handle FCSV format
    header <- readLines(landmark_path, n = 3)
    lmks <- read.table(landmark_path, sep = ",", skip = 3, header = TRUE)

    # Apply mirroring
    lmks[, "x"] <- lmks[, "x"] * x
    lmks[, "y"] <- lmks[, "y"] * y
    lmks[, "z"] <- lmks[, "z"] * z

    # Overwrite file with header and new landmark coordinates
    writeLines(header, con = landmark_path)
    write.table(lmks, file = landmark_path, append = TRUE, sep = ",",
                row.names = FALSE, col.names = TRUE)

  } else if (file_ext == "json") {
    # Handle JSON format
    data <- read_json(landmark_path)  # no simplifyVector!

    for (i in seq_along(data$markups[[1]]$controlPoints)) {
      coords <- as.numeric(data$markups[[1]]$controlPoints[[i]]$position)
      data$markups[[1]]$controlPoints[[i]]$position <- c(coords[1] * x,
                                                         coords[2] * y,
                                                         coords[3] * z)
    }

    write_json(data, landmark_path, auto_unbox = TRUE, pretty = TRUE)

  } else {
    stop("Unsupported file type. Only .fcsv and .json formats are supported.")
  }
}

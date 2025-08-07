test_that("bone_scan_check works with mesh input", {
  skip_if_not_installed("curl")
  skip_if_not(curl::has_internet(), "No internet connection")

  # STL mesh file
  url_mesh <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
  mesh_file <- tempfile(fileext = ".stl")
  mesh_success <- tryCatch({
    download.file(url_mesh, mesh_file, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Failed to download mesh: ", e$message)
    FALSE
  })

  skip_if_not(mesh_success, "Mesh download failed")

  surface_mesh <- tryCatch({
    import_mesh(mesh_file)
  }, error = function(e) {
    skip(paste("import_mesh failed:", e$message))
  })

  # CT scan file
  url_scan <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
  scan_file <- tempfile(fileext = ".nii.gz")
  scan_success <- tryCatch({
    download.file(url_scan, scan_file, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Failed to download scan: ", e$message)
    FALSE
  })

  skip_if_not(scan_success, "Scan download failed")

  nifti <- tryCatch({
    import_scan(scan_file)
  }, error = function(e) {
    skip(paste("import_scan failed:", e$message))
  })

  df <- bone_scan_check(surface_mesh, nifti, return_limits = TRUE)

  expect_s3_class(df, "data.frame")
  expect_true(all(c("Axis", "Mesh_Min", "Mesh_Max", "Scan_Min", "Scan_Max") %in% colnames(df)))
})

test_that("bone_scan_check works with matrix input", {
  skip_if_not_installed("curl")
  skip_if_not(curl::has_internet(), "No internet connection")

  # STL mesh file
  url_mesh <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
  mesh_file <- tempfile(fileext = ".stl")
  mesh_success <- tryCatch({
    download.file(url_mesh, mesh_file, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Failed to download mesh: ", e$message)
    FALSE
  })

  skip_if_not(mesh_success, "Mesh download failed")

  surface_mesh <- tryCatch({
    import_mesh(mesh_file)
  }, error = function(e) {
    skip(paste("import_mesh failed:", e$message))
  })

  # CT scan file
  url_scan <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
  scan_file <- tempfile(fileext = ".nii.gz")
  scan_success <- tryCatch({
    download.file(url_scan, scan_file, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Failed to download scan: ", e$message)
    FALSE
  })

  skip_if_not(scan_success, "Scan download failed")

  nifti <- tryCatch({
    import_scan(scan_file)
  }, error = function(e) {
    skip(paste("import_scan failed:", e$message))
  })

  coords <- t(surface_mesh$vb)[, 1:3]
  df <- bone_scan_check(coords, nifti, return_limits = TRUE)

  expect_s3_class(df, "data.frame")
  expect_true(all(c("Axis", "Mesh_Min", "Mesh_Max", "Scan_Min", "Scan_Max") %in% colnames(df)))
})

test_that("bone_scan_check works with data.frame input", {
  skip_if_not_installed("curl")
  skip_if_not(curl::has_internet(), "No internet connection")

  # STL mesh file
  url_mesh <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
  mesh_file <- tempfile(fileext = ".stl")
  mesh_success <- tryCatch({
    download.file(url_mesh, mesh_file, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Failed to download mesh: ", e$message)
    FALSE
  })

  skip_if_not(mesh_success, "Mesh download failed")

  surface_mesh <- tryCatch({
    import_mesh(mesh_file)
  }, error = function(e) {
    skip(paste("import_mesh failed:", e$message))
  })

  coords_df <- as.data.frame(t(surface_mesh$vb)[, 1:3])

  # CT scan file
  url_scan <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
  scan_file <- tempfile(fileext = ".nii.gz")
  scan_success <- tryCatch({
    download.file(url_scan, scan_file, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Failed to download scan: ", e$message)
    FALSE
  })

  skip_if_not(scan_success, "Scan download failed")

  nifti <- tryCatch({
    import_scan(scan_file)
  }, error = function(e) {
    skip(paste("import_scan failed:", e$message))
  })

  result <- bone_scan_check(coords_df, nifti)
  expect_null(result)  # returns invisible(NULL) when return_limits = FALSE
})

test_that("bone_scan_check errors when mesh extends outside scan volume", {
  skip_if_not_installed("curl")
  skip_if_not(curl::has_internet(), "No internet connection")

  # STL mesh file
  url_mesh <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
  mesh_file <- tempfile(fileext = ".stl")
  mesh_success <- tryCatch({
    download.file(url_mesh, mesh_file, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Failed to download mesh: ", e$message)
    FALSE
  })

  skip_if_not(mesh_success, "Mesh download failed")

  surface_mesh <- tryCatch({
    import_mesh(mesh_file)
  }, error = function(e) {
    skip(paste("import_mesh failed:", e$message))
  })

  coords_df <- as.data.frame(t(surface_mesh$vb)[, 1:3])

  # CT scan file
  url_scan <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
  scan_file <- tempfile(fileext = ".nii.gz")
  scan_success <- tryCatch({
    download.file(url_scan, scan_file, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Failed to download scan: ", e$message)
    FALSE
  })

  skip_if_not(scan_success, "Scan download failed")

  nifti <- tryCatch({
    import_scan(scan_file)
  }, error = function(e) {
    skip(paste("import_scan failed:", e$message))
  })

  mesh_outside <- surface_mesh
  mesh_outside$vb[1, ] <- mesh_outside$vb[1, ] + 10000

  expect_error(
    bone_scan_check(mesh_outside, nifti, return_limits = TRUE),
    "Mesh not within scan volume."
  )
})

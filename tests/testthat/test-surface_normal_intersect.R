landmark_path <- system.file("extdata", "test_femur.mrk.json", package = "BoneDensityMapping")
landmarks <- import_lmks(landmark_path)

test_that("surface_normal_intersect returns expected output", {
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

  # Generate surface points
  mapped_coords <- surface_points_template(surface_mesh, landmarks, no_surface_sliders = 10)

  # Run the function
  mat_peak <- surface_normal_intersect(
    surface_mesh, mapped_coords, normal_dist = 3.0,
    nifti
  )

  # Check output type and length
  expect_length(mat_peak, nrow(mapped_coords))

  # Check for NAs (optional depending on what you expect)
  expect_false(any(is.na(mat_peak)))
})

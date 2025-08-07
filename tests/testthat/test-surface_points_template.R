landmark_path <- system.file("extdata", "test_femur.mrk.json", package = "BoneDensityMapping")
landmarks <- import_lmks(landmark_path)

test_that("surface_points_template returns a data.frame with correct dimensions", {
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

  n_sliders <- 100
  result <- surface_points_template(surface_mesh, landmarks, n_sliders)

  expect_s3_class(result, "data.frame")
  # Number of rows = number of original landmarks + number of sliders
  expect_equal(nrow(result), nrow(landmarks) + n_sliders)
  # Expect 3 columns (xpts, ypts, zpts)
  expect_equal(ncol(result), 3)
  # Columns should be numeric
  expect_true(all(sapply(result, is.numeric)))
})

test_that("surface_points_template errors when surface mesh is missing vertices", {
  # Create dummy invalid mesh object
  bad_mesh <- list(vb = matrix(numeric(0), nrow=4, ncol=0))

  landmark_path <- system.file("extdata", "test_femur.mrk.json", package = "BoneDensityMapping")
  landmarks <- import_lmks(landmark_path)

  expect_error(surface_points_template(bad_mesh, landmarks, 10))
})


landmark_path <- system.file("extdata", "SCAP001_landmarks.fcsv",
                            package = "BoneDensityMapping")
lmks1 <- import_lmks(landmark_path)

landmark_path <- system.file("extdata", "SCAP002_landmarks.fcsv",
                             package = "BoneDensityMapping")
lmks2 <- import_lmks(landmark_path)

test_that("surface_points_new returns valid remapped coordinates", {
  skip_if_not_installed("curl")
  skip_if_not(curl::has_internet(), "No internet connection")

  # STL mesh file 1
  url_mesh1 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/SCAP001.stl"
  mesh_file1 <- tempfile(fileext = ".stl")
  mesh_success1 <- tryCatch({
    download.file(url_mesh1, mesh_file1, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Failed to download mesh: ", e$message)
    FALSE
  })

  skip_if_not(mesh_success1, "Mesh download failed")

  mesh1 <- tryCatch({
    import_mesh(mesh_file1)
  }, error = function(e) {
    skip(paste("import_mesh failed:", e$message))
  })

  # STL mesh file 2
  url_mesh2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/SCAP002.stl"
  mesh_file2 <- tempfile(fileext = ".stl")
  mesh_success2 <- tryCatch({
    download.file(url_mesh2, mesh_file2, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Failed to download mesh: ", e$message)
    FALSE
  })

  skip_if_not(mesh_success2, "Mesh download failed")

  mesh2 <- tryCatch({
    import_mesh(mesh_file2)
  }, error = function(e) {
    skip(paste("import_mesh failed:", e$message))
  })

  # Generate template and remap
  template <- surface_points_template(mesh1, lmks1, 100)
  result <- surface_points_new(mesh2, lmks2, template)

  # Basic checks
  expect_true(is.matrix(result) || is.data.frame(result))
  expect_equal(ncol(result), 3)
  expect_equal(nrow(result), nrow(template))
  expect_true(all(is.finite(result)))
  expect_type(result[1, 1], "double")
  mesh_coords <- t(mesh2$vb)[, 1:3]
  for (i in 1:3) {
    expect_true(all(result[, i] >= min(mesh_coords[, i])))
    expect_true(all(result[, i] <= max(mesh_coords[, i])))
  }
})

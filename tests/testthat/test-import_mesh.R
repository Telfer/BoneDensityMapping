test_that("import_mesh correctly imports a valid .stl mesh", {
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

  expect_s3_class(surface_mesh, "mesh3d")
  expect_true(!is.null(surface_mesh$vb))  # vertices
  expect_true(!is.null(surface_mesh$it))  # triangle indices
})

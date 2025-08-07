test_that("fill_bone_points returns correct point matrix", {
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

  # Run fill_bone_points with spacing 10
  points <- fill_bone_points(surface_mesh, spacing = 10)

  # Test 1: output is matrix
  expect_true(is.matrix(points))

  # Test 2: points are within bounding box of mesh vertices (plus small margin)
  verts <- t(surface_mesh$vb)[, 1:3]
  margin <- 0.01
  expect_true(all(points[,1] >= min(verts[,1]) - margin))
  expect_true(all(points[,1] <= max(verts[,1]) + margin))
  expect_true(all(points[,2] >= min(verts[,2]) - margin))
  expect_true(all(points[,2] <= max(verts[,2]) + margin))
  expect_true(all(points[,3] >= min(verts[,3]) - margin))
  expect_true(all(points[,3] <= max(verts[,3]) + margin))

  # Test 3: some points returned
  expect_gt(nrow(points), 0)

  # Test 4: number of points returned is not greater than total points in grid
  x <- seq(min(verts[,1]), max(verts[,1]), by = 10)
  y <- seq(min(verts[,2]), max(verts[,2]), by = 10)
  z <- seq(min(verts[,3]), max(verts[,3]), by = 10)
  expect_lte(nrow(points), length(x) * length(y) * length(z))
})

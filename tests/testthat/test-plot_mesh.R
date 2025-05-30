test_that("plot_mesh runs without errors", {
  surface_mesh_path <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")
  surface_mesh <- import_mesh(surface_mesh_path)

  # Minimal test with just required arguments
  expect_silent(
    plot_mesh(surface_mesh, title = "Test Mesh")
  )
})

test_that("plot_mesh works with density colors and user matrix", {
  surface_mesh_path <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")
  surface_mesh <- import_mesh(surface_mesh_path)

  # Create dummy color and user matrix
  n_vertices <- ncol(surface_mesh$vb)
  density_color <- rep("red", n_vertices)

  userMat <- diag(4)  # identity matrix as a dummy orientation

  expect_silent(
    plot_mesh(surface_mesh, title = "Test Colored Mesh", density_color = density_color, userMat = userMat)
  )
})

test_that("mesh_template_match returns correct nearest neighbor indices", {
  # Create a fake mesh with a few vertices
  # Note: surface_mesh$vb is a 4xN matrix; last row is typically 1s (homogeneous coordinates)
  vb <- rbind(
    c(0, 1, 0, 0),
    c(0, 0, 1, 0),
    c(0, 0, 0, 1),
    c(1, 1, 1, 1)
  )
  surface_mesh <- list(vb = vb)

  # Create template points (each row = a template point)
  template_points <- matrix(c(
    0, 0, 0,   # Closest to vertex 1
    1, 0, 0,   # Closest to vertex 2
    0, 1, 0,   # Closest to vertex 3
    0, 0, 1    # Closest to vertex 4
  ), ncol = 3, byrow = TRUE)

  # Run the function
  matched <- mesh_template_match(surface_mesh, template_points)

  # Expect a vector of length = number of mesh vertices (4)
  expect_type(matched, "integer")
  expect_length(matched, 4)

  # Each vertex matches exactly to the corresponding template point
  expect_equal(matched, 1:4)
})

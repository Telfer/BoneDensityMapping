test_that("import_mesh correctly imports a valid .stl mesh", {
  path <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")
  mesh <- import_mesh(path)

  expect_true(inherits(mesh, "mesh3d"))
  expect_true(!is.null(mesh$vb))  # vertices
  expect_true(!is.null(mesh$it))  # triangle indices
})

test_that("import_mesh errors on unsupported file format", {
  fake_file <- tempfile(fileext = ".txt")
  writeLines("not a mesh", fake_file)

  expect_error(import_mesh(fake_file))
})


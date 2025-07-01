test_that("surface_points_new returns valid remapped coordinates", {
  # Load mesh and landmark data
  mesh1_path <- system.file("extdata", "SCAP001.stl", package = "BoneDensityMapping")
  lmk1_path <- system.file("extdata", "SCAP001_landmarks.fcsv", package = "BoneDensityMapping")
  mesh2_path <- system.file("extdata", "SCAP002.stl", package = "BoneDensityMapping")
  lmk2_path <- system.file("extdata", "SCAP002_landmarks.fcsv", package = "BoneDensityMapping")

  mesh1 <- import_mesh(mesh1_path)
  lmks1 <- import_lmks(lmk1_path)
  mesh2 <- import_mesh(mesh2_path)
  lmks2 <- import_lmks(lmk2_path)

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

library(rgl)

test_that("surface_points_new looks reasonable", {
  skip_if_not(interactive())  # prevent in CI

  surface_mesh_path <- system.file("extdata", "SCAP001.stl",
                                   package = "BoneDensityMapping")
  scap_001_mesh <- import_mesh(surface_mesh_path)
  landmark_path <- system.file("extdata", "SCAP001_landmarks.fcsv",
                                package = "BoneDensityMapping")
  scap_001_lmk <- import_lmks(landmark_path)
  template_coords <- surface_points_template(scap_001_mesh, scap_001_lmk, 1000)

  surface_mesh_path <- system.file("extdata", "SCAP002.stl",
                                    package = "BoneDensityMapping")
  scap_002_mesh <- import_mesh(surface_mesh_path)
  landmark_path <- system.file("extdata", "SCAP002_landmarks.fcsv",
                               package = "BoneDensityMapping")
  scap_002_lmk <- import_lmks(landmark_path)
  scap_002_remapped <- surface_points_new(scap_002_mesh, scap_002_lmk, template_coords)

  open3d() +
    shade3d(scap_001_mesh, color = "gray", alpha = 0.2) +
    points3d(scap_001_lmk[, c("x", "y", "z")], col = "blue", size = 2) +
    points3d(template_coords, col = "green", size = 0.5) +

    shade3d(scap_002_mesh, color = "gray", alpha = 0.2) +
    points3d(scap_002_lmk[, c("x", "y", "z")], col = "red", size = 2) +
    points3d(scap_002_remapped, col = "purple", size = 0.5)

  expect(TRUE)
})

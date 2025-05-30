test_that("bone_scan_check works with mesh input", {
  scan_path <- system.file("extdata", "test_CT_hip.nii", package = "BoneDensityMapping")
  mesh_path <- system.file("extdata", "test_CT_femur.stl",package = "BoneDensityMapping")

  nifti <- import_scan(scan_path)
  mesh <- import_mesh(mesh_path)

  # Expect no error for valid mesh inside scan
  expect_output(bone_scan_check(mesh, nifti), "Axis")
})

test_that("bone_scan_check works with matrix input", {
  scan_path <- system.file("extdata", "test_CT_hip.nii", package = "BoneDensityMapping")
  mesh_path <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")

  nifti <- import_scan(scan_path)
  mesh <- import_mesh(mesh_path)

  coords <- t(mesh$vb)[, 1:3]  # extract vertex coordinates
  expect_output(bone_scan_check(coords, nifti), "Axis")
})

test_that("bone_scan_check works with data.frame input", {
  scan_path <- system.file("extdata", "test_CT_hip.nii", package = "BoneDensityMapping")
  mesh_path <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")

  nifti <- import_scan(scan_path)
  mesh <- import_mesh(mesh_path)

  coords_df <- as.data.frame(t(mesh$vb)[, 1:3])
  expect_output(bone_scan_check(coords_df, nifti), "Axis")
})

test_that("bone_scan_check errors on unsupported surface_mesh types", {
  scan_path <- system.file("extdata", "test_CT_hip.nii", package = "BoneDensityMapping")
  nifti <- import_scan(scan_path)

  expect_error(bone_scan_check(list(1, 2, 3), nifti),
               "surface_mesh must be a mesh3d object or a matrix of vertex coordinates")
  expect_error(bone_scan_check(NULL, nifti),
               "surface_mesh must be a mesh3d object or a matrix of vertex coordinates")
})

test_that("bone_scan_check errors when mesh extends outside scan volume", {
  scan_path <- system.file("extdata", "test_CT_hip.nii", package = "BoneDensityMapping")
  mesh_path <- system.file("extdata", "test_CT_femur.stl",package = "BoneDensityMapping")

  nifti <- import_scan(scan_path)
  mesh <- import_mesh(mesh_path)

  # Artificially shift mesh vertices outside volume
  mesh_outside <- mesh
  mesh_outside$vb[1, ] <- mesh_outside$vb[1, ] + 10000  # Shift X coords far out of range

  expect_error(
    bone_scan_check(mesh_outside, nifti),
    "Mesh not within scan volume"
  )
})


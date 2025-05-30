test_that("voxel_point_intersect returns correct length vector", {
  mesh_path <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")
  mesh <- import_mesh(mesh_path)

  nifti_path <- system.file("extdata", "test_CT_hip.nii", package = "BoneDensityMapping")
  nifti <- import_scan(nifti_path)

  # Use a small number of points
  vertices <- t(mesh$vb)[1:5, 1:3]

  # Run function
  mat_peak <- voxel_point_intersect(vertices, nifti, betaCT = 1.0, sigmaCT = 2.0)

  # Check output
  expect_type(mat_peak, "double")
  expect_length(mat_peak, nrow(vertices))
  expect_false(any(is.na(mat_peak)))
})

test_that("voxel_point_intersect runs with check_in_vol = TRUE", {
  mesh_path <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")
  mesh <- import_mesh(mesh_path)

  nifti_path <- system.file("extdata", "test_CT_hip.nii", package = "BoneDensityMapping")
  nifti <- import_scan(nifti_path)

  vertices <- t(mesh$vb)[1:5, 1:3]

  expect_no_error({
    voxel_point_intersect(vertices, nifti, betaCT = 1.0, sigmaCT = 2.0, check_in_vol = TRUE)
  })
})


test_that("surface_normal_intersect returns expected output", {
  # Load test data
  surface_mesh_path <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")
  surface_mesh <- import_mesh(surface_mesh_path)
  landmark_path <- system.file("extdata", "test_femur.mrk.json", package = "BoneDensityMapping")
  landmarks <- import_lmks(landmark_path)
  nifti_path <- system.file("extdata", "test_CT_hip.nii", package = "BoneDensityMapping")
  nifti <- import_scan(nifti_path)

  # Generate surface points
  mapped_coords <- surface_points_template(surface_mesh, landmarks, no_surface_sliders = 10)

  # Run the function
  mat_peak <- surface_normal_intersect(
    surface_mesh, mapped_coords, normal_dist = 3.0,
    nifti, betaCT = 1.0, sigmaCT = 1.0
  )

  # Check output type and length
  expect_type(mat_peak, "double")
  expect_length(mat_peak, nrow(mapped_coords))

  # Check for NAs (optional depending on what you expect)
  expect_false(any(is.na(mat_peak)))
})

test_that("surface_normal_intersect errors if mesh is out of scan bounds", {
  # Load test data
  surface_mesh_path <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")
  surface_mesh <- import_mesh(surface_mesh_path)
  landmark_path <- system.file("extdata", "test_femur.mrk.json", package = "BoneDensityMapping")
  landmarks <- import_lmks(landmark_path)
  nifti_path <- system.file("extdata", "test_CT_hip.nii", package = "BoneDensityMapping")
  nifti <- import_scan(nifti_path)

  mapped_coords <- surface_points_template(surface_mesh, landmarks, no_surface_sliders = 10)

  # Intentionally shift mapped_coords far outside volume to trigger error
  mapped_coords_shifted <- mapped_coords
  mapped_coords_shifted[, c("xpts", "ypts", "zpts")] <- mapped_coords_shifted[, c("xpts", "ypts", "zpts")] + 1000

  expect_error(
    surface_normal_intersect(
      surface_mesh, mapped_coords_shifted, normal_dist = 3.0,
      nifti, betaCT = 1.0, sigmaCT = 1.0
    ),
    "Mesh not within scan volume"
  )
})

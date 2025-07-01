test_that("landmark_check passes when all landmarks are within threshold", {
  mesh_path <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")
  lmk_path <- system.file("extdata", "test_femur.fcsv", package = "BoneDensityMapping")

  mesh <- import_mesh(mesh_path)
  landmarks <- import_lmks(lmk_path)

  expect_message(
    landmark_check(mesh, landmarks, threshold = 2.0),
    "Landmarks are on bone surface")
})


test_that("landmark_check warns when landmarks exceed threshold", {
  mesh_path <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")
  lmk_path <- system.file("extdata", "test_femur.fcsv", package = "BoneDensityMapping")

  mesh <- import_mesh(mesh_path)
  landmarks <- import_lmks(lmk_path)

  expect_message(
    landmark_check(mesh, landmarks, threshold = 0.5),
    "Landmarks not on bone surface"
  )
})

test_that("landmark_check errors with malformed landmarks", {
  mesh_path <- system.file("extdata", "test_CT_femur.stl", package = "BoneDensityMapping")
  mesh <- import_mesh(mesh_path)

  # Missing coordinate columns
  bad_landmarks <- data.frame(lmk_id = c("A", "B"))

  expect_error(landmark_check(mesh, bad_landmarks))
})


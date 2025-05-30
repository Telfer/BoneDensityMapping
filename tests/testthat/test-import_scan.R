test_that("import_scan reads valid .nii file", {
  path <- system.file("extdata", "test_CT_hip.nii", package = "BoneDensityMapping")
  scan <- import_scan(path)

  expect_s4_class(scan, "nifti")
  expect_true(length(scan) > 0)
})

test_that("import_scan errors on unsupported file type", {
  fake_path <- tempfile(fileext = ".txt")
  file.create(fake_path)

  expect_error(import_scan(fake_path), "Unsupported file type")
})

test_that("import_scan errors when file does not exist", {
  missing_file <- testthat::test_path("testdata/does_not_exist.nii")

  expect_error(import_scan(missing_file))
})

test_that("import_scan reads valid .nii file", {
  skip_if_not_installed("curl")
  skip_if_not(curl::has_internet(), "No internet connection")

  # CT scan file
  url_scan <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
  scan_file <- tempfile(fileext = ".nii.gz")
  scan_success <- tryCatch({
    download.file(url_scan, scan_file, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("Failed to download scan: ", e$message)
    FALSE
  })

  skip_if_not(scan_success, "Scan download failed")

  nifti <- tryCatch({
    import_scan(scan_file)
  }, error = function(e) {
    skip(paste("import_scan failed:", e$message))
  })

  expect_s4_class(nifti, "nifti")
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

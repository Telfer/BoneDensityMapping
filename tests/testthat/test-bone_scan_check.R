url <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.2/test_CT_femur.stl"
bone_filepath <- tempfile(fileext = ".stl")
download.file(url, bone_filepath, mode = "wb")
surface_mesh <- import_mesh(bone_filepath)

url2 <- "https://github.com/Telfer/BoneDensityMapping/releases/download/v1.0.1/test_CT_hip.nii.gz"
scan_filepath <- tempfile(fileext = ".nii.gz")
download.file(url2, scan_filepath, mode = "wb")
nifti <- import_scan(scan_filepath)

test_that("bone_scan_check works with mesh input", {

  # Expect no error for valid mesh inside scan
  expect_output(bone_scan_check(surface_mesh, nifti), "Axis")
})

test_that("bone_scan_check works with matrix input", {
  coords <- t(surface_mesh$vb)[, 1:3]  # extract vertex coordinates
  expect_output(bone_scan_check(coords, nifti), "Axis")
})

test_that("bone_scan_check works with data.frame input", {


  coords_df <- as.data.frame(t(surface_mesh$vb)[, 1:3])
  expect_output(bone_scan_check(coords_df, nifti), "Axis")
})

test_that("bone_scan_check errors on unsupported surface_mesh types", {

  expect_error(bone_scan_check(list(1, 2, 3), nifti),
               "surface_mesh must be a mesh3d object or a matrix of vertex coordinates")
  expect_error(bone_scan_check(NULL, nifti),
               "surface_mesh must be a mesh3d object or a matrix of vertex coordinates")
})

test_that("bone_scan_check errors when mesh extends outside scan volume", {

  # Artificially shift mesh vertices outside volume
  mesh_outside <- surface_mesh
  mesh_outside$vb[1, ] <- mesh_outside$vb[1, ] + 10000  # Shift X coords far out of range

  expect_error(
    bone_scan_check(mesh_outside, nifti),
    "Mesh not within scan volume"
  )
})


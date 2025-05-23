% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoneDensityFunctions.R
\name{fill_bone_points}
\alias{fill_bone_points}
\title{fill bone}
\usage{
fill_bone_points(surface_mesh, spacing)
}
\arguments{
\item{surface_mesh}{Mesh object}

\item{spacing}{Numeric}
}
\value{
Matrix with internal point coordinates
}
\description{
fill bone
}
\examples{
surface_mesh_path <- system.file("extdata", "test_CT_femur.stl",
                                 package = "BoneDensityMapping")
surface_mesh <- import_mesh(surface_mesh_path)
internal_fill <- fill_bone_points(surface_mesh, 10)
}

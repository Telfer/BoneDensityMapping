% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoneDensityFunctions.R
\name{surface_points_new}
\alias{surface_points_new}
\title{New surface points from template}
\usage{
surface_points_new(surface_mesh, landmarks, template)
}
\arguments{
\item{surface_mesh}{List. Mesh data imported via ply_import function}

\item{landmarks}{Data frame. Contains 3D coords of landmarks}

\item{template}{Data frame. 3D coords of remapped surface points}
}
\value{
Data frame. 3D coords of remapped surface points
}
\description{
New surface points from template
}
\author{
Scott Telfer \email{scott.telfer@gmail.com} Adpated from geomorph
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoneDensityFunctions.R
\name{rm_local_sig}
\alias{rm_local_sig}
\title{local significance}
\usage{
rm_local_sig(vertices, sig_vals, sig_level = 0.05, dist, n_local = 1)
}
\arguments{
\item{vertices}{Matrix}

\item{sig_vals}{Numeric vector}

\item{sig_level}{Numeric. Default 0.05}

\item{dist}{Numeric. Distance to check for vertices}

\item{n_local}{Numeric. Number of local significant values needed}
}
\value{
Numeric vector
}
\description{
local significance
}
\author{
Scott Telfer \email{scott.telfer@gmail.com}
}

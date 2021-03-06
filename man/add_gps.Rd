% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_prepare_data.R
\name{add_gps}
\alias{add_gps}
\title{Add generalized PS}
\usage{
add_gps(data, formula, family, subset, ps_prefix)
}
\arguments{
\item{data}{data_frame}

\item{formula}{formula for PS model}

\item{family}{\code{multinomial(parallel = FALSE)} should be used.}

\item{subset}{subset expression if estimating PS only for the subset}

\item{ps_prefix}{string used as the prefix. Use to distinguish the re-estimated PS from the original PS.}
}
\value{
data_frame containing three additional columns for PS.
}
\description{
.. content for details ..
}

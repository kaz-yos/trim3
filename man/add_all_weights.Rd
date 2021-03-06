% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_prepare_data.R
\name{add_all_weights}
\alias{add_all_weights}
\title{Calculate one type of PS weight given PS}
\usage{
add_all_weights(data, A_name, levels, ps_prefix1 = "ps1_",
  ps_prefix2 = "ps2_")
}
\arguments{
\item{data}{data_frame}

\item{levels}{Character vector holding three elements corresponding to levels in \code{A}}

\item{ps_prefix1}{Prefix for the PS estimated in the entire cohort.}

\item{ps_prefix2}{Prefix for the PS estimated in the trimmed cohort.}
}
\value{
data_frame containing IPTW, MW, and OW estimated in the entire cohort as well as the trimmed cohort.
}
\description{
.. content for details ..
}

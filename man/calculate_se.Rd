% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference.R
\name{calculate_se}
\alias{calculate_se}
\title{Calculate Standard Errors Using Multiplier Bootstrap}
\usage{
calculate_se(
  inf_matrix,
  cluster_id = NULL,
  biter = 2000,
  pl = TRUE,
  n_cores = 8,
  alp = 0.05
)
}
\arguments{
\item{inf_matrix}{NxK matrix with N indiv observations and K ATTs}

\item{cluster_id}{Vector nx1 of cluster IDs.}

\item{biter}{Number of bootstrap iterations}

\item{pl}{Run in parallel}

\item{n_cores}{Number of cores to use in parallel}

\item{alp}{Test size, defaults to 0.05.}
}
\description{
Calculate standard errors using influence function
}

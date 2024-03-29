% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference.R
\name{run_multiplier_bootstrap}
\alias{run_multiplier_bootstrap}
\title{Run Multiplier Bootstrap}
\usage{
run_multiplier_bootstrap(inf.func, biters, pl = FALSE, cores = 1, ...)
}
\arguments{
\item{inf.func}{Influence function matrix. NxK where N is N indiv and K is # of ATTs}

\item{biters}{Bootstrap iterations. N.B. plural use of arg here but biter elsewhere. Nice.}

\item{pl}{Process in parallel}

\item{cores}{Number of cores to use if parallel processing}

\item{...}{Unused argument, just added to keep calls to run_nested_multiplier_bootstrap the same

TAKEN DIRECTLY FROM BCALLAWAY11/DID}
}
\description{
Multiplier Bootstrap
}

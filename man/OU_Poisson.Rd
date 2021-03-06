% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonGaussianOU.R
\name{OU_Poisson}
\alias{OU_Poisson}
\title{OU_Poisson}
\usage{
OU_Poisson(Tmax, Nsim, x0, alpha, mu)
}
\arguments{
\item{Tmax}{time at which we want to sample from;}

\item{Nsim}{number of simulated values from Y(Tmax)}

\item{x0}{non-negative starting condition}

\item{alpha}{positive exponential decay rate}

\item{mu}{positive intensity of the underlying Poisson process}
}
\description{
Simulation of Nsim values from Y(Tmax), with Y(t) being non-Gaussian OU-Poisson (see Tamborrino and Lansky, Physica D, 418 (2021) 132845, TL2021),
i.e. an OU process driven by a Poisson process with jump 1 and intensity mu>0
}

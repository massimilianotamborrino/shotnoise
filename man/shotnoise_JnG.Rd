% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shotnoise_sim.r
\name{shotnoise_JnG}
\alias{shotnoise_JnG}
\title{shotnoise_JnG}
\usage{
shotnoise_JnG(Tmax, Nsim, x0, alpha, mu, tilde.sigma2, lambda)
}
\arguments{
\item{Tmax}{time at which we want to sample from;}

\item{Nsim}{number of simulated values from Xn(Tmax)}

\item{x0}{non-negative starting condition}

\item{alpha}{positive exponential decay rate}

\item{mu}{positive quantity entering into the scale and rate parameters of the Gamma jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
This corresponds to the limit first infinitesimal moment, eq. (23)}

\item{tilde.sigma2}{positive quantity entering into the scale and rate parameters of the Gamma jump amplitude of Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
This corresponds to the limit second infinitesimal moment, eq. (24), see also Table 1}

\item{lambda}{positive intensity of the underlying Poisson process}
}
\value{
Nsim values from Xn(Tmax) of the shot noise process with jumps Jn being Gamma distributed
}
\description{
Simulation of Nsim values from the shot noise process Xn(t) with Gamma distributed jumps Jn at time Tmax
}

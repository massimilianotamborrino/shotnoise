% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonGaussianOU.R
\name{OU_IG}
\alias{OU_IG}
\title{OU_IG}
\usage{
OU_IG(Tmax, Ntime, Nsim, x0, alpha, mu, tilde.sigma2)
}
\arguments{
\item{Tmax}{time at which we want to sample from;}

\item{Ntime}{number of equidistant points to simulate in [0,Tmax]}

\item{Nsim}{number of simulated values from Y(t)}

\item{x0}{non-negative starting condition}

\item{alpha}{positive exponential decay rate}

\item{mu}{positive quantity entering into the mean and shape parameters of the IG jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
This corresponds to the limit first infinitesimal moment, eq. (23)}

\item{tilde.sigma2}{positive quantity entering into the shape parameter  of the IG jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
This corresponds to the limit second infinitesimal moment, eq. (24), see also Table 1}
}
\description{
Simulation of Nsim values from Y(Tmax), with Y(t) being non-Gaussian OU-Inverse Gaussian (OU-IG) process (13) of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021),
using the algorithm in Y. Qu, A. Dassios, H. Zhao, J. Appl. Probab. (2021) [QDZ2021].
dY(t) = -alpha Y(t)dt + 1*dZ(t), with Z(t) being an Inverse-Gaussian process with
Lévy measure nu(dx)=exp(-c^2 x/2)/sqrt(2*pi*x^3)dx, with c=sqrt(mu/tilde.sigma^2), rho=s=(from [QDZ2021])= sqrt(mu^3/tilde.sigma^2)
 (see Section 3.2 and Table 1 of TL2021)
}

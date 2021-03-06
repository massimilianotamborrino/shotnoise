% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonGaussianOU.R
\name{OU_Gamma}
\alias{OU_Gamma}
\title{OU_Gamma}
\usage{
OU_Gamma(Tmax, Ntime, Nsim, x0, alpha, mu, tilde.sigma2, scenario)
}
\arguments{
\item{Tmax}{time at which we want to sample from;}

\item{Ntime}{total number of time discretisation points of the OU-Gamma;}

\item{Nsim}{number of simulated values;}

\item{x0}{starting condition of the processes;}

\item{alpha}{positive exponential decay rate of the shot noise process}

\item{mu}{positive quantity entering into the scale and rate parameters of the Gamma jump amplitude Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
This corresponds to the limit first infinitesimal moment, eq. (23)}

\item{tilde.sigma2}{positive quantity entering into the scale and rate parameters of the Gamma jump amplitude of Jn, see Table 1 of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021).
This corresponds to the limit second infinitesimal moment, eq. (23), see also Table 1}

\item{scenario}{1 if beta=0.5 (OU-Gamma obtained as limiting process of a shot noise with Chi-square distributed jumps) or 2 if beta=mu/tilde.sigma^2 (OU-Gamma obtained as limiting process of a shot noise with Gamma distributed jumps)}
}
\value{
Nsim simulated values of the OU-Gamma process Y at time Tmax
}
\description{
Simulation of Nsim values from Y(Tmax), with Y(t) being non-Gaussian OU-Gamma process (13) of Tamborrino and Lansky, Physica D, 418 (2021) 132845 (TL2021),
using the algorithm in Y. Qu, A. Dassios, H. Zhao, J. Oper. Res. Soc. (2019)
dY(t) = -alpha Y(t)dt + 1*dZ(t), with Z(t) being a gamma process with
Lévy measure nu(dx)=tilde.alpha x^(-1)+exp(-beta x)dx,
with (tilde.alpha,beta)=
(mu/2,1/2) if the OU-Gamma process has been obtained as limit of shot noise with chi-square distributed jumps,
or (mu^2/tilde.sigma^2, mu/tilde.sigma^2) if the OU-Gamma process has been obtained as limit of shot noise with Gamma distributed jumps;
 (see Section 3.2 and Table 1 of TL2021)
}

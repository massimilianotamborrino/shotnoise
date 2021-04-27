# shotnoise
R package for the exact simulation of non-negative shot noise processes (with Bernoulli, Poisson, Gamma and inverse-Gaussian distributed jumps) and non-Gaussian OU processes, also known as Lévy-driven OU processes (in particular, OU-Poisson, OU-Gamma and OU-inverse Gaussian (OU-IG) processes) from the paper

[1] M. Tamborrino, P. Lansky, Shot noise, weak convergence and diffusion approximation, Physica D, 418, 132845, 2021, https://www.sciencedirect.com/science/article/abs/pii/S0167278921000038

The code was written by Massimiliano Tamborrino (firstname dot secondname at warwick.ac.uk) in Rcpp. Yan Qu provided us with the Matlab codes for the exact simulation of the OU-Gamma and OU-IG processes based on the papers

[2] Y. Qu, A. Dassios, H. Zhao, Exact simulation of Gamma-driven Ornstein–
Uhlenbeck processes with finite and infinite activity jumps, J. Oper. Res.
Soc., 471--484, 2019, http://dx.doi.org/10.1080/01605682.2019.1657368.

[3] Y. Qu, A. Dassios, H. Zhao, Exact simulation of Ornstein–Uhlenbeck
tempered stable processes, J. Appl. Probab., 2021 (in press). Preprint at http://eprints.lse.ac.uk/106267/

Those Matlab algorithms have been then rewritten in Rcpp and included in the package.

# What can you find in the package
Here we provide the code for
1) the exact simulation of N values from a non-negative shot noise process $X$ at time $t$, solution of the stochastic differential equation (2) from [1]

dX(t)= -alpha X(t) dt + J dN(t),

where N(t) is a Poisson process with rate lambda>0 and J is the distribution of the jump amplitude. In the package, we consider J to be a Bernoulli distribution (shotnoise_JnBer), Poisson (shotnoise_JnP), chi-square (shotnoise_JnP with parameter scenario equal to 1), Gamma (shotnoise_JnG, with scenario equal to 2), inverse-Gaussian ( shotnoise_JnIG). The chosen parameters and the generating algorithm is described in [1].
2) the exact simulation of a trajectory of the shot noise process with Gamma and IG jump amplitudes (shotnoise_JnG_trajectory and shotnoise_JnIG_trajectory, respectively).
3) the exact simulation of $N$ values from an OU-Poisson process at time $t$, solution of (13) of [1]
\begin{equation}\label{13} 
dY(t)= -\delta Y(t) dt + \rho dZ(t)\end{equation},
with $\delta=\alpha$ and $Z(t)$ being a Poisson process with intensity $\mu>0$, see Section 3 of [1]. R code: OU_Poisson
4) the exact simulation of $N$ values from an OU-Gamma process at time $t$, solution of (13) of [1] 
$$ dY(t)= -\delta Y(t) dt + \rho dZ(t)$$,
with $\delta=\alpha$ and $Z(t)$ being a gamma process with Lévy measure $\nu(dx)=\tilde \alpha x^{-1}e^{-\beta x}dx$. See Section 3.2 and Table 1 for the chosen parameters, and [2] for the description of the simulation algorithm. R code: OU_Gamma
5) the exact simulation of $N$ values from an OU-IG process at time $t$, solution of (13) of [1] 
$$ dY(t)= -\delta Y(t) dt + \rho dZ(t)$$,
with $\delta=\alpha$ and $Z(t)$ being an Inverse Gaussian process with Lévy measure $\nu(dx)=e^{-c^2 x/2}dx$, with $c^2=\mu/\tilde\sigma^2$. See Section 3.2 and Table 1 for the chosen parameters, and [3] for the description of the simulation algorithm.  R code: OU_IG
6) the simulation of a trajectory of the OU-Gamma process (OU_Gamma_trajectory) and the OU-IG process (OU_IG_trajectory)
7) the computation of integrated absolute errors (IAEs) between the shot noise process and the corresponding Gaussian OU or non-Gaussian OU process.

# How to install the package
The simplest way is to do it via devtools, using devtools::install_github("massimilianotamborrino/shotnoise")

# How does the package work
Use code_Figure1.r to generate two trajectories of the OU-Gamma and OU-IG processes, reproducing Fig. 1 of [1].

Use code_Figure2.r to compare the estimated probability density functions of the shot noise process, the Gaussian OU and the non-Gaussian OU, reproducing Fig. 2 of [1]. 

All parameters are chosen according to Table 1 of [1].

# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

OU_IG_V_generator2 <- function(w) {
    .Call('_shotnoise_OU_IG_V_generator2', PACKAGE = 'shotnoise', w)
}

function_IG_generator <- function(n, mu, lambda) {
    .Call('_shotnoise_function_IG_generator', PACKAGE = 'shotnoise', n, mu, lambda)
}

OU_IG_Rcpp_ <- function(T_terminal, N_time, N_path, X_0, alpha, mu, sigma2) {
    .Call('_shotnoise_OU_IG_Rcpp_', PACKAGE = 'shotnoise', T_terminal, N_time, N_path, X_0, alpha, mu, sigma2)
}

OU_IG_trajectory_Rcpp_ <- function(T_terminal, N_time, X_0, alpha, mu, sigma2) {
    .Call('_shotnoise_OU_IG_trajectory_Rcpp_', PACKAGE = 'shotnoise', T_terminal, N_time, X_0, alpha, mu, sigma2)
}

stl_sortb <- function(x) {
    .Call('_shotnoise_stl_sortb', PACKAGE = 'shotnoise', x)
}

OU_Poisson_Rcpp_ <- function(t0, Nsim, x0, alpha, mu) {
    .Call('_shotnoise_OU_Poisson_Rcpp_', PACKAGE = 'shotnoise', t0, Nsim, x0, alpha, mu)
}

OU_Gamma_Rcpp_ <- function(T_terminal, N_time, N_path, Y_0, delta, mu, sigma2, scenario) {
    .Call('_shotnoise_OU_Gamma_Rcpp_', PACKAGE = 'shotnoise', T_terminal, N_time, N_path, Y_0, delta, mu, sigma2, scenario)
}

OU_Gamma_trajectory_Rcpp_ <- function(T_terminal, N_time, Y_0, delta, mu, sigma2, scenario) {
    .Call('_shotnoise_OU_Gamma_trajectory_Rcpp_', PACKAGE = 'shotnoise', T_terminal, N_time, Y_0, delta, mu, sigma2, scenario)
}

stl_sort <- function(x) {
    .Call('_shotnoise_stl_sort', PACKAGE = 'shotnoise', x)
}

cumsum_sug <- function(x) {
    .Call('_shotnoise_cumsum_sug', PACKAGE = 'shotnoise', x)
}

shotnoise_JnP_Rcpp_ <- function(t0, Nsim, x0, alpha, mu, lambda) {
    .Call('_shotnoise_shotnoise_JnP_Rcpp_', PACKAGE = 'shotnoise', t0, Nsim, x0, alpha, mu, lambda)
}

shotnoise_JnBer_Rcpp_ <- function(t0, Nsim, x0, alpha, mu, lambda) {
    .Call('_shotnoise_shotnoise_JnBer_Rcpp_', PACKAGE = 'shotnoise', t0, Nsim, x0, alpha, mu, lambda)
}

shotnoise_JnG_Rcpp_ <- function(t0, Nsim, x0, alpha, mu, sigma2, lambda) {
    .Call('_shotnoise_shotnoise_JnG_Rcpp_', PACKAGE = 'shotnoise', t0, Nsim, x0, alpha, mu, sigma2, lambda)
}

shotnoise_JnG_trajectory_Rcpp_ <- function(t0, N_time, x0, alpha, mu, sigma2, lambda) {
    .Call('_shotnoise_shotnoise_JnG_trajectory_Rcpp_', PACKAGE = 'shotnoise', t0, N_time, x0, alpha, mu, sigma2, lambda)
}

shotnoise_JnChiSq_Rcpp_ <- function(t0, Nsim, x0, alpha, mu, lambda) {
    .Call('_shotnoise_shotnoise_JnChiSq_Rcpp_', PACKAGE = 'shotnoise', t0, Nsim, x0, alpha, mu, lambda)
}

shotnoise_JnIG_Rcpp_ <- function(t0, Nsim, x0, alpha, mu, sigma2, lambda) {
    .Call('_shotnoise_shotnoise_JnIG_Rcpp_', PACKAGE = 'shotnoise', t0, Nsim, x0, alpha, mu, sigma2, lambda)
}

shotnoise_JnIG_trajectory_Rcpp_ <- function(t0, N_time, x0, alpha, mu, sigma2, lambda) {
    .Call('_shotnoise_shotnoise_JnIG_trajectory_Rcpp_', PACKAGE = 'shotnoise', t0, N_time, x0, alpha, mu, sigma2, lambda)
}


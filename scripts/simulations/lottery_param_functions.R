# References for the biological choices below: scripts/simulations/references.bib

#' Draw positions along a life-history axis
#'
#' Each species gets a value \eqn{z \in [0, 1]} on a perennial--annual axis,
#' where \eqn{z = 0} is a long-lived, low-fecundity, short-dispersing perennial
#' and \eqn{z = 1} is a short-lived, high-fecundity, far-dispersing annual.
#' All other species-level parameters are derived from \eqn{z} by
#' \code{life_history_params()}.
#'
#' @param nsp Number of species
#' @param shape1,shape2 Beta shape parameters used when \code{even = FALSE}.
#' \code{shape1 < shape2} skews the community toward perennials.
#' @param even If \code{TRUE}, species are spaced evenly along the axis
#' (\code{ppoints(nsp)}); otherwise \eqn{z} is drawn from a Beta distribution.
#'
#' @return Vector of length \code{nsp} on the unit interval, sorted from most
#' perennial to most annual.
#' @export
#'
#' @examples
#' draw_life_history(10)
#' draw_life_history(10, shape1 = 1, shape2 = 3, even = FALSE)
#'
draw_life_history <- function(nsp, shape1 = 1, shape2 = 1, even = TRUE){
  if(even){
    z <- ppoints(nsp)
  } else{
    z <- rbeta(nsp, shape1 = shape1, shape2 = shape2)
  }
  return(sort(z))
}








#' Species parameters from the life-history axis
#'
#' Maps the life-history axis \eqn{z} onto death probability, fecundity,
#' dispersal, environmental sensitivity, environmental optima, and the strength
#' of intraspecific competition.
#'
#' Fecundity is tied to death rate (\eqn{\lambda_i \propto \delta_i}) as an
#' equalizing trade-off: in a lottery model a rare species invades a resident
#' when \eqn{\lambda_i / \delta_i > \lambda_r / \delta_r}, so holding
#' \eqn{\lambda / \delta} roughly constant keeps fitness differences small and
#' lets stabilizing mechanisms (intra > inter competition, local crowding, the
#' storage effect, immigration) maintain many species. Seeds that stay in the
#' parent's (occupied) cell are wasted, so fecundity is divided by
#' \eqn{1 - P(d = 0)} to equalize the seeds that actually reach other cells.
#' \code{fitness_sd} sets how far species depart from equal fitness.
#'
#' @param z Vector of life-history positions (0 = perennial, 1 = annual)
#' @param sigma_env Standard deviation of the innovations in the AR(1) environment
#' @param phi Autocorrelation of the AR(1) environment
#' @param d_max Max dispersal distance (as passed to \code{seed_rain_array_fast()})
#' @param p_death Death probability for \eqn{z = 0} and \eqn{z = 1}
#' @param disp_rate Exponential dispersal rate for \eqn{z = 0} and \eqn{z = 1}.
#' Higher rates keep seeds near the parent (clumping).
#' @param f0 Max fecundity (seeds leaving the parent cell) of a species with
#' death probability 1
#' @param fitness_sd SD of log-scale noise around the fecundity--death trade-off
#' @param sigma_env_sens Environmental sensitivity (\code{sigma} in
#' \code{fecundity_ll()}) for \eqn{z = 0} and \eqn{z = 1}
#' @param alpha Mean strength of intraspecific competition
#' @param alpha_sd SD of log-scale noise in intraspecific competition
#'
#' @return Named list of per-species vectors: \code{z}, \code{Pr_death},
#' \code{lambda_max}, \code{disp_rate}, \code{sigma}, \code{optims}, \code{alpha}.
#' @export
#'
#' @examples
#' life_history_params(draw_life_history(5), sigma_env = 0.1, phi = 0.4, d_max = 5)
#'
life_history_params <- function(z, sigma_env, phi, d_max,
                                p_death = c(0.03, 1),
                                disp_rate = c(2, 0.2),
                                f0 = 150,
                                fitness_sd = 0.1,
                                sigma_env_sens = c(10, 40),
                                alpha = 4,
                                alpha_sd = 0.05){
  nsp <- length(z)

  # long-lived block at low z, strict annual at z = 1
  Pr_death <- p_death[1] + (p_death[2] - p_death[1]) * z^2

  # log-linear interpolation of dispersal rate (short for perennials)
  rate <- exp((1 - z) * log(disp_rate[1]) + z * log(disp_rate[2]))

  # equalizing trade-off: seeds leaving the parent cell proportional to death rate
  p_stay <- sapply(rate, function(x){P_d_exp(d = 0, d_max = d_max, lambda = x)})
  lambda_max <- f0 * Pr_death / (1 - p_stay) * exp(rnorm(nsp, 0, fitness_sd))

  # annuals track the environment more closely
  sigma <- sigma_env_sens[1] + diff(sigma_env_sens) * z

  # optima spread over the stationary env distribution, independent of z
  env_sd <- sigma_env / sqrt(1 - phi^2)
  optims <- sample(qnorm(ppoints(nsp), sd = env_sd))

  list(
    z = z,
    Pr_death = Pr_death,
    lambda_max = lambda_max,
    disp_rate = rate,
    sigma = sigma,
    optims = optims,
    alpha = alpha * exp(rnorm(nsp, 0, alpha_sd))
  )
}








#' Competition matrix along a life-history gradient
#'
#' Builds an \eqn{S \times S} matrix where \eqn{A_{ij}} is the competitive
#' effect of species \eqn{j} on the fecundity of species \eqn{i}:
#' \deqn{A_{ij} = \alpha_i \rho_{max} \exp\left(-\frac{(z_i - z_j)^2}{2 w^2}\right)
#' \exp\left(h \frac{z_i - z_j}{2}\right) e^{\epsilon_{ij}}, \quad A_{ii} = \alpha_i.}
#' The first factor is limiting similarity (similar life histories compete more);
#' the second makes competition asymmetric, with perennials (low \eqn{z})
#' suppressing annuals more than the reverse.
#'
#' @param z Vector of life-history positions
#' @param alpha Vector of intraspecific competition strengths
#' @param rho_max Interspecific competition between identical life histories,
#' relative to intraspecific competition
#' @param overlap_width Width \eqn{w} of the limiting-similarity kernel on \eqn{z}
#' @param hierarchy Strength \eqn{h} of the perennial-over-annual asymmetry
#' (0 = symmetric)
#' @param noise_sd SD of log-scale noise \eqn{\epsilon_{ij}}
#'
#' @return \eqn{S \times S} matrix of competition coefficients
#' @export
#'
#' @examples
#' z <- draw_life_history(5)
#' comp_matrix_gradient(z, alpha = rep(4, 5))
#'
comp_matrix_gradient <- function(z, alpha, rho_max = 0.5, overlap_width = 0.3,
                                 hierarchy = 1, noise_sd = 0.1){
  nsp <- length(z)
  dz <- outer(z, z, "-")

  overlap <- exp(-dz^2 / (2 * overlap_width^2))
  asym <- exp(hierarchy * dz / 2)
  noise <- matrix(exp(rnorm(nsp^2, 0, noise_sd)), nrow = nsp)

  mat <- alpha * rho_max * overlap * asym * noise
  diag(mat) <- alpha

  if(any(mat > diag(mat) & row(mat) != col(mat))){
    warning("some interspecific coefficients exceed intraspecific ones")
  }

  return(mat)
}








#' Add immigrant seeds from a regional species pool
#'
#' Adds a uniform trickle of seeds to every cell, split among species by
#' \code{meta_abund}. The total immigrant rain per cell is \code{m} times the
#' mean total local seed rain per cell, so \code{m} is a relative rate. This
#' keeps species from going extinct across the whole landscape and guarantees
#' every cell receives some seeds.
#'
#' @param seed_rain \eqn{M \times J \times S} seed rain array
#' @param m Immigration rate relative to mean local seed rain (0 = closed)
#' @param meta_abund Relative abundance of each species in the regional pool
#' (defaults to equal)
#'
#' @return Seed rain array with immigrant seeds added
#' @export
#'
#' @examples
#' sr <- array(rgamma(48, 2), dim = c(4, 4, 3))
#' add_immigration(sr, m = 0.01)
#'
add_immigration <- function(seed_rain, m, meta_abund = NULL){
  S <- dim(seed_rain)[3]
  if(is.null(meta_abund)){
    meta_abund <- rep(1 / S, S)
  }
  meta_abund <- meta_abund / sum(meta_abund)

  # floor at 1 seed per cell so the trickle never vanishes
  per_cell <- m * max(1, sum(seed_rain) / prod(dim(seed_rain)[1:2]))

  for(s in 1:S){
    seed_rain[, , s] <- seed_rain[, , s] + per_cell * meta_abund[s]
  }
  return(seed_rain)
}








#' Coexistence and clumping diagnostics for a simulated landscape
#'
#' @param ts_all List of lattices (one per time step)
#' @param z Vector of life-history positions
#' @param burn_in Number of initial time steps to ignore when averaging cover
#'
#' @return Data frame with one row per species: \code{z}, mean cover after
#' burn-in, final cover, whether the species persists at the final step, and a
#' clumping index (mean fraction of conspecifics among the 8 nearest neighbors
#' divided by landscape cover; > 1 means clumped).
#' @export
#'
#' @examples
#'
landscape_diagnostics <- function(ts_all, z, burn_in = 0){
  nsp <- length(z)
  keep <- ts_all[(burn_in + 1):length(ts_all)]
  X_T <- ts_all[[length(ts_all)]]

  cover_t <- sapply(keep, function(X){tabulate(X, nbins = nsp) / length(X)})
  final_cover <- tabulate(X_T, nbins = nsp) / length(X_T)

  # fraction of each cell's neighbors that are conspecific, at r = 1
  nbrs <- kernel_count_fast(X_T, r = 1, sp_list = 1:nsp)
  n_nbrs <- apply(nbrs, c(1, 2), sum)
  clumping <- sapply(1:nsp, function(s){
    occ <- X_T == s
    if(!any(occ)) return(NA)
    mean(nbrs[, , s][occ] / n_nbrs[occ]) / final_cover[s]
  })

  out <- data.frame(
    species = 1:nsp,
    z = z,
    mean_cover = rowMeans(cover_t),
    final_cover = final_cover,
    persists = final_cover > 0,
    clumping = clumping
  )
  attr(out, "n_persist") <- sum(out$persists)
  return(out)
}

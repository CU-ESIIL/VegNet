# document setup
library(here)
library(ggplot2)

source(here("scripts/simulations/terrestrial_lottery_sim_functions.R"))
source(here("scripts/simulations/lottery_param_functions.R"))

# setup parameters

  # lattice size along one margin
  cells <- 60

  # neighborhood size
  nbrhood_radius <- 3

  # max dispersal distance
  d_max <- 5

  # time steps (the first burn_in steps are discarded after the run)
  steps <- 300
  burn_in <- 100

  # immigration from the regional pool, relative to mean local seed rain
  m <- 0.01

  # environment
  sigma_env <- 0.1    # variation in the environment
  phi <- 0.4          # autocorrelation across time steps

  # species placed along a perennial (z = 0) to annual (z = 1) axis
  nsp <- 15
  z <- draw_life_history(nsp)

  # death, fecundity, dispersal, env response, and self-regulation from z
  pars <- life_history_params(z, sigma_env = sigma_env, phi = phi, d_max = d_max)
  Pr_death <- pars$Pr_death
  lambda_max <- pars$lambda_max
  disp_rate <- pars$disp_rate
  sp_optims <- pars$optims

  # competition matrix: limiting similarity + perennial-over-annual asymmetry
  A_mat <- comp_matrix_gradient(z, alpha = pars$alpha)

  # initialize lattice
  X <- matrix(
    data = sample(1:nsp, cells * cells, replace = T),
    nrow = cells, ncol = cells
  )

  # count the neighbors for each cell given the radius and
  #  lattice size
  n_neighbors <- count_neighbors(M = cells, J = cells, r = nbrhood_radius)

  # track evolution of community
  ts_all <- vector(mode = "list", length = steps)
  ts_all[[1]] <- X

  # track environment
  env <- vector(mode = "double", length = steps)
  env[1] <- rnorm(1, sd = sigma_env/sqrt(1 - phi^2))


  # step through the evolution of the community
  for(t in 1:(steps-1)){

    # count neighbors of each species around each cell
    nbrs_t <- kernel_count_fast(ts_all[[t]], r = nbrhood_radius, sp_list = 1:nsp)

    # get fecundity matrix by applying the fecundidty_ll() command to each row of X
    Fecundity <- t(sapply(
      1:cells,
      FUN = function(x, lattice, alpha, nbrhood, n){
        fecundity_ll(
          foc_sp = lattice[x, ],
          lambda_max = lambda_max,
          optims = sp_optims,
          env_t = env[t],
          alpha = alpha,
          nbrhood = nbrhood[x, ],
          n = n[x, , ],
          sigma = pars$sigma
        )
      },
      lattice = ts_all[[t]],
      alpha = A_mat,
      nbrhood = n_neighbors,
      n = nbrs_t
    ))

    # calculate seed rain into each cell by each species
    sr_t <- seed_rain_array_fast(
      F_mat = Fecundity,
      X = ts_all[[t]],
      d_max = d_max,
      rate = disp_rate,
      nsp = nsp
    )

    # add immigrant seeds from the regional pool
    sr_t <- add_immigration(sr_t, m = m)

    # kill off and replace some adults with seedlings
    ts_all[[t+1]] <- die_replace(
      ts_all[[t]],
      prob_death = Pr_death,
      seed_rain = sr_t
    )

    # change environment for next iteration
    env[t + 1] <- phi * env[t] + rnorm(1, sd = sigma_env)

  }

# coexistence and clumping along the life-history axis (after burn-in)
  diagnostics <- landscape_diagnostics(ts_all, z = z, burn_in = burn_in)
  diagnostics
  attr(diagnostics, "n_persist")

# drop the burn-in so ts_all[[t]] and env[t] stay aligned downstream
  ts_all <- ts_all[(burn_in + 1):steps]
  env <- env[(burn_in + 1):steps]
  steps <- length(ts_all)

# get percent cover dataframe
  cover_df <- data.frame(
    t = rep(1:steps, nsp),
    species = rep(1:nsp, each = steps)
  )

# fill in the cover values
  cover <- vector(mode = "double")
  for(i in 1:nsp){
    cover <- c(
      cover,
      sapply(ts_all, function(x){sum(x == i)/(nrow(x) * ncol(x))})
    )
  }
  cover_df$cover <- cover
  cover_df$z <- z[cover_df$species]
  cover_df$species <- as.factor(cover_df$species)

  # color by life history: 0 = perennial, 1 = annual
  ggplot(data = cover_df, aes(x = t, y = cover, group = species, color = z))+
    geom_line()+
    scale_color_viridis_c()+
    theme_classic()






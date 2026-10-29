# document setup
library(here)
library(ggplot2)

source(here("scripts/simulations/terrestrial_lottery_sim_functions.R"))

# setup parameters

  # lattice size along one margin
  cells <- 50

  # neighborhood size
  nbrhood_radius <- 3

  # time steps
  steps <- 50

  # species
  S <- 5 # number of strong competitors
  s <- 0 # number weak competitors
  nsp <- S + s

  # vector of intra-specific competitive effects
#  alpha <- rgamma(nsp, shape = 3, rate = 1)
  alpha <- c(
    sort(runif(S+1, min = 3.8, max = 4)),
    runif(s, min = 3.8, max = 4)
  )

  # competition matrix
  # A_mat <- comp_matrix(n_sp = nsp, alpha = alpha, rho = c(0.01, 0.4), num_ngs = S, num_regs = 0)
  A_mat <- comp_matrix2(n_sp = nsp, alpha = alpha, rho = 0, num_ngs = S)


  # max per-capita fecundity
  # lambda <- c(
  #   rgamma(S + 1, shape = 40, rate = 1),
  #   rgamma(s, shape = 10, rate = 1)
  # )
  lambda_max <- c(
    runif(5, min = 45, max = 50),
    runif(45, min = 10, max = 15)
  )
  sp_optims <- runif(nsp, min = -0.2, max = 0.2)

  # probability an individual dies in a given time point for each species
  # Pr_death <- c(
  #   rep(1, 1 + S),
  #   runif(s, min = 0.1, max = 0.5)
  # )
  Pr_death <- c(
    rep(1, 5),
    runif(45, min = 0.05, max = 0.2)
  )

  # dispersal rates for each species
  disp_rate <- c(
    runif(S),
    runif(s, min = 0, max = 3)
  )

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
  sigma_env <- 0.1    # variation in the environment
  phi <- 0.4          # autocorrelation across time steps
  env <- vector(mode = "double", length = steps)
  env[1] <- rnorm(1, sd = sigma_env/sqrt(1 - phi^2))


  # step through the evolution of the community
  for(t in 1:(steps-1)){

    # count neighbors of each species around each cell
    nbrs_t <- kernel_count(ts_all[[t]], r = nbrhood_radius, sp_list = 1:nsp)

    # get fecundity matrix by applying the fecundidty_ll() command to each row of X
    Fecundity <- t(sapply(
      1:cells,
      FUN = function(x, lattice, lambda_t, alpha, nbrhood, n){
        fecundity_ll(
          foc_sp = lattice[x, ],
          lambda_max = lambda_max,
          optims = sp_optims,
          env_t = env[t],
          alpha = alpha,
          nbrhood = nbrhood[x, ],
          n = n[x, , ]
        )
      },
      lattice = ts_all[[t]],
      lambda_t = lambda,
      alpha = A_mat,
      nbrhood = n_neighbors,
      n = nbrs_t
    ))

    # calculate seed rain into each cell by each species
    sr_t <- seed_rain_array(
      F_mat = Fecundity,
      X = ts_all[[t]],
      d_max = 5,
      rate = disp_rate,
      nsp = nsp
    )

    # kill off and replace some adults with seedlings
    ts_all[[t+1]] <- die_replace(
      ts_all[[t]],
      prob_death = Pr_death,
      seed_rain = sr_t
    )

    # change environment for next iteration
    env[t + 1] <- phi * env[t] + rnorm(1, sd = sigma_env)

  }

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
  cover_df$species <- as.factor(cover_df$species)

  ggplot(data = cover_df, aes(x = t, y = cover, color = species))+
    geom_line()+
    theme_classic()+
    theme(legend.position = "none")

# are the strongly competing species negatively correlated?
  cover_trunc <- cover_df[cover_df$t >= 100, ]
  cor.test(
    cover_trunc$cover[cover_trunc$species == 2],
    cover_trunc$cover[cover_trunc$species == 3]
  )

# how many species coexist?
  sum(cover_df[cover_df$t == steps, ]$cover > 0)






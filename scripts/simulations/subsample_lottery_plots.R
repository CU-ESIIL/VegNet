# Subsample a simulated lottery-model landscape into plots and structure
#  the plot-level cover data in long, lagged format:
#  plot_id | time | species | env | cover_tm1 | cover_t

# document setup
library(here)
library(ggplot2)

# randomly choose n_plots non-overlapping square sub-grids of side
#  plot_size from a cells x cells lattice
sample_plot_corners <- function(cells, plot_size, n_plots){

  # tile the lattice into candidate blocks
  starts <- seq(1, cells - plot_size + 1, by = plot_size)
  blocks <- expand.grid(row_start = starts, col_start = starts)

  if(n_plots > nrow(blocks)){
    stop(
      "n_plots (", n_plots, ") exceeds the ", nrow(blocks),
      " non-overlapping plots that fit in the lattice"
    )
  }

  corners <- blocks[sample(nrow(blocks), n_plots), ]
  data.frame(
    plot_id = 1:n_plots,
    row_start = corners$row_start,
    col_start = corners$col_start
  )
}

# proportion cover of each species in each plot at each time step
#  (zeros kept explicitly)
plot_cover <- function(ts_all, corners, plot_size, nsp){

  out <- vector(mode = "list", length = nrow(corners))
  for(p in 1:nrow(corners)){

    rows <- corners$row_start[p] + 0:(plot_size - 1)
    cols <- corners$col_start[p] + 0:(plot_size - 1)

    # time x species matrix of cover for plot p
    cover_p <- t(sapply(
      ts_all,
      function(X){tabulate(X[rows, cols], nbins = nsp) / plot_size^2}
    ))

    out[[p]] <- data.frame(
      plot_id = corners$plot_id[p],
      time = rep(1:length(ts_all), nsp),
      species = rep(1:nsp, each = length(ts_all)),
      cover = as.vector(cover_p)
    )
  }

  do.call(rbind, out)
}

# pair cover at t (t = 2...T) with cover at t - 1 (t = 1...T-1)
lag_cover <- function(cover_long, env){

  lagged <- cover_long
  lagged$time <- lagged$time + 1
  names(lagged)[names(lagged) == "cover"] <- "cover_tm1"

  # inner join drops t = 1 (no lag) and t = T + 1 (no response)
  out <- merge(cover_long, lagged, by = c("plot_id", "species", "time"))
  names(out)[names(out) == "cover"] <- "cover_t"

  # environment that drives the t-1 -> t transition, matching
  #  X_env[p, t - 1] in scripts/stan/ricker_ss_hurdle_lognormal.stan
  out$env <- env[out$time - 1]

  out$species <- as.factor(out$species)
  out <- out[order(out$plot_id, out$species, out$time), ]
  rownames(out) <- NULL
  out[, c("plot_id", "time", "species", "env", "cover_tm1", "cover_t")]
}

# subsample the simulated landscape

  # number of plots and side length of each (in cells)
  n_plots <- 30
  plot_size <- 5

  # run the simulation to get ts_all, env, nsp, cells, and steps
  #  NOTE: raise `cells` in lottery_model_sims.R (e.g., 100) so the
  #  landscape fits n_plots non-overlapping plots
  source(here("scripts/simulations/lottery_model_sims.R"))

  corners <- sample_plot_corners(cells = cells, plot_size = plot_size, n_plots = n_plots)
  cover_long <- plot_cover(ts_all, corners = corners, plot_size = plot_size, nsp = nsp)
  plot_df <- lag_cover(cover_long, env = env)

  head(plot_df)

# sanity check: cover at t against cover at t - 1
  ggplot(data = plot_df, aes(x = cover_tm1, y = cover_t, color = env))+
    geom_point(alpha = 0.4)+
    facet_wrap(~species)+
    theme_classic()

# save locally (keep large simulated data out of GitHub)
  # saveRDS(plot_df, "~/vegnet_sims/lottery_plot_df.rds")

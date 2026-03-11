simdata_filename <- function(cfg, data_dir) {
  sprintf(
    "%s%s_simdata_n%g_sd%.1f_MC%g_seed%g.RData",
    data_dir,
    cfg$experiment,
    cfg$n,
    cfg$sd,
    cfg$N_MC,
    cfg$seed
  )
}

load_simulation_data <- function(cfg, data_dir) {
  file <- simdata_filename(cfg, data_dir)
  
  if (!file.exists(file)) {
    stop("Simulation data file not found:\n", file)
  }
  
  message("Loading data from:\n", file)
  load(file, envir = .GlobalEnv)
  data = data[, 1:cfg$n]
  #invisible(file)
  return(data)
}
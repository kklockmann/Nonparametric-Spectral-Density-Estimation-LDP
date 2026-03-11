source("plotting_utils.R")
source("ldp_estimation.R")
source("data_utils.R")

#Before running 
#1. simulate data, which are loaded here
#2. create two folders: FIGURES and RESULTS
data_dir = "../DATA/"
fig_dir = "../FIGURES/"
res_dir = "../RESULTS/"

#Configuration
cfg <- list(
  experiment = "EX1",
  #switch between "EX1", "EX2", "EX3", "EX4"
  
  #index of covariance coefficient to be estimated
  j = 2,  #(0: variance, 1: simga_1, ...)
  
  #data parameters
  n = 1000,
  sd = 1.2,
  
  #Monte Carlo
  N_MC = 300,
  seed = 42,
  
  #theory parameters
  delta = 0.1,
  
  #to reduce theoretical choice of tau.tilde:
  #theoretical choice: 1;
  #other choices, e.g., 16, 160
  division_const_tilde = 160
)

###########################################
# Computation for covariance estimation
###########################################
#initialize
set.seed(cfg$seed)

# -----------------------
#Load data
# -----------------------
data = load_simulation_data(cfg, data_dir)

#true covariance coefficient
acf_true = acf.true[cfg$j + 1] #indexing: acf.true[1]=sigma_0,...

#grid of privacy levels
alpha_grid =  exp(seq(log(0.01), log(100000), length.out = 100))

#theoretical truncation parameters
tau_NI_theo <- sqrt(56 * log(cfg$n)^(1 + cfg$delta))
tau_SI_theo <- sqrt(8 * log(cfg$n)^(1 + cfg$delta))


#list of truncation parameters to be used
#Note: labels in legend of the plots have to changed manually below,
#      in case other division factors, or more truncation parameters
#      are use than the default
if (cfg$j == 0) {
  tau_NI <- c(tau_NI_theo, tau_NI_theo / 6.5, tau_NI_theo / 20)
  tau_SI <- c(tau_SI_theo, tau_SI_theo / 1.5, tau_SI_theo / 7)
} else{
  tau_NI <- c(tau_NI_theo, tau_NI_theo / 6.5, tau_NI_theo / 20) #5,3,2,1)
  tau_SI <- c(tau_SI_theo, tau_SI_theo / 2.5, tau_SI_theo / 7) #5,3,2,1)
}
tau_tilde_SI <- 16 * log(cfg$n)^(1 + cfg$delta) * tau_SI^2 / cfg$division_const_tilde

tau_grid_acf = cbind(tau_NI, tau_SI, tau_tilde_SI)



# ---------------------------
# No privacy (np) baseline
# ---------------------------
mse_np <- numeric(cfg$N_MC)
for (i in seq_len(cfg$N_MC)) {
  X <- data[i, ]
  sigma_hat <- mean(X[1:(cfg$n - cfg$j)] * X[(cfg$j + 1):cfg$n])
  mse_np[i] <- (sigma_hat - acf_true)^2
}
MSE_np <- mean(mse_np)


# -----------------------
# Private estimators
# -----------------------
MSE_SI <- MSE_NI <- matrix(0, length(alpha_grid), length(tau_NI))

for (a in 1:length(alpha_grid)) {
  alpha = alpha_grid[a]
  
  for (t in 1:nrow(tau_grid_acf)) {
    tau <- tau_grid_acf[t, 1]
    tau_si <- tau_grid_acf[t, 2]
    tau_tilde <- tau_grid_acf[t, 3]
    
    mse_SI <- mse_NI <- numeric(cfg$N_MC)
    
    for (i in seq_len(cfg$N_MC)) {
      X <- data[i, ]
      
      # Sequential interactive privatization and estimation
      Z_SI <- priv.data.SI.acf(X, cfg$j, alpha, tau_si, tau_tilde)
      mse_SI[i] <- (sigma.hat.SI(Z_SI) - acf_true)^2
      
      # Non-interactive privatization and estimation
      Z_NI <- priv.data.NI(X, alpha, tau)
      mse_NI[i] <- (sigma.hat.NI(Z_NI, cfg$j, alpha, tau) - acf_true)^2
    }
    MSE_SI[a, t] <- mean(mse_SI)
    MSE_NI[a, t] <- mean(mse_NI)
  }
}

results.SI <- results.NI <- list(
  target = "acf",
  alpha_grid = alpha_grid,
  tau_grid = tau_grid_acf,
  MSE_p = MSE_NI,
  MSE_np = MSE_np,
  acf_true = acf_true,
  cfg = cfg,
  typ = "NI"
)
results.SI$MSE_p =  MSE_SI
results.SI$typ = "SI"

#up to which alpha has SI-based estimator smaller MSE?
#max(alpha_grid[results.SI$MSE_p[,1]<= results.NI$MSE_p[,1]])#theoretical tau
#max(alpha_grid[results.SI$MSE_p[,2]<= results.NI$MSE_p[,2]])#chosen tau


###############################################################
#Plots for covariance estimation: log(MSE) versus log(alpha)
###############################################################
#NI-plot: variance and covariance
#manually adapt labels in case division factors for tau have been changed above
labels_NI = c(
  expression(tau == tau[theory]),
  bquote(tau == tau[theory] / 6.5),
  bquote(tau == tau[theory] / 20),
  "no privacy",
  expression(-4 * log(alpha))
)

plot_mse_vs_alpha(
  results.NI,
  fig_dir,
  log_scale_y = TRUE,
  log_scale_x = TRUE,
  theo_rate = "NI",
  labels = labels_NI,
  ylims = c(-5,30)
)

#SI-plot: variance
#manually adapt labels in case division factors for tau have been changed above
if (cfg$j == 0) {
  labels_SI = c(
    expression(tau == tau[theory]),
    bquote(tau == tau[theory] / 1.5),
    bquote(tau == tau[theory] / 7),
    "no privacy",
    expression(-2 * log(alpha))
  )
  
  #ylims = vector with lowest and highest value on y-axis,
  #        in order to manually align axis of different plots
  plot_mse_vs_alpha(
    results.SI,
    fig_dir,
    log_scale_y = TRUE,
    log_scale_x = TRUE,
    theo_rate = "SI",
    labels = labels_SI,
    ylims = c(-5,30)
  )
}

#SI-plot: covariance estimation for j not equal to 0
#manually adapt labels in case division factors for tau have been changed above
if (cfg$j != 0) {
  labels_SI = c(
    bquote(
      tau == tau[theory] * "," ~  ~ tilde(tau) == tilde(tau)[theory] / .(cfg$division_const_tilde)
    ),
    bquote(
      tau == tau[theory] / 2.5 * "," ~  ~ tilde(tau) == tilde(tau)[theory] / .(cfg$division_const_tilde)
    ),
    bquote(
      tau == tau[theory] / 7 * "," ~  ~ tilde(tau) == tilde(tau)[theory] / .(cfg$division_const_tilde)
    ),
    "no privacy",
    expression(-2 * log(alpha))
  )
  
  #ylims = vector with lowest and highest value on y-axis,
  #        in order to manually align axis of different plots
  plot_mse_vs_alpha(
    results.SI,
    fig_dir,
    log_scale_y = TRUE,
    log_scale_x = TRUE,
    theo_rate = "SI",
    labels = labels_SI,
    ylims = c(-5,30)
  )
}

#save results
results_file =  
  sprintf(
    "%s%s_acf_n%d_j%d_MC%d_seed%d.RDATA",
    res_dir,
    cfg$experiment,
    cfg$n,
    cfg$j,
    cfg$N_MC,
    cfg$seed
  )

save(cfg, results.SI, results.NI, file = results_file)

# #######################################################################
# #for j!=0 - Plot in paper
# #with two different division factors for tau.tilde
# #######################################################################
# #1. run the above code with selected j
# #   and cfg$division_const_tilde = 160 (example of first choice)
# #2. set
# results.SI_ti1 = results.SI
# #3. run the above code with same j
# #   and cfg$division_const_tilde = 1 (example of second choice)
# #   without running the line results.SI_ti1 = results.SI
# #4. execute
# results.SI_ti1$MSE_p = cbind(results.SI$MSE_p[, 1], results.SI_ti1$MSE_p[, 2:3])
# results.SI_ti1$tau_grid = rbind(results.SI$tau_grid[1, ], results.SI_ti1$tau_grid[2:3, ])
# #5. Plot (eventually adapt division factors manually if default was changed)
# labels_SI = c(
#   expression(tau == tau[theory] * "," ~  ~ tilde(tau) == tilde(tau)[theory]),
#   bquote(tau == tau[theory] / 2.5 * "," ~  ~ tilde(tau) == tilde(tau)[theory] / 160),
#   bquote(tau == tau[theory] / 7 * "," ~  ~ tilde(tau) == tilde(tau)[theory] / 160),
#   "no privacy",
#   expression(-2 * log(alpha))
# )
# 
# #ylims = vector with lowest and highest value on y-axis,
# #        in order to manually align axis of different plots
# plot_mse_vs_alpha(
#   results.SI_ti1,
#   fig_dir,
#   log_scale_y = TRUE,
#   log_scale_x = TRUE,
#   theo_rate = "SI",
#   labels = labels_SI,
#   ylims = c(-10,30)
# )
#
#save(cfg,results.SI_ti1,results.SI,results.NI,file=results_file)

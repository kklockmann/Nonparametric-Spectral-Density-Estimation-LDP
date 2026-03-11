source("plotting_utils.R")
source("ldp_estimation.R")
source("data_utils.R")

#Before running 
#1. simulate data, which are loaded here
#2. create two folders: FIGURES and RESULTS
data_dir = "../DATA/"
fig_dir = "../FIGURES/"
res_dir = "../RESULTS/"

#to speed up the program, use shorter alpha_grid in computations section

#Configuration
cfg <- list(
  experiment = "EX1",
  #switch between "EX1", "EX2", "EX3", "EX4"
  
  #frequency of sdf to be estimated, 
  omega = pi / 5, #value in [0,pi]
  
  #data parameters
  n = 1000,
  sd = 1.2,
  
  #Monte Carlo
  N_MC = 300,
  seed = 42,
  
  #theory parameters
  delta = 0.001,
  s = 3,  #smoothness of sdf
  
  #to reduce theoretical choice of tau.tilde:
  #theoretical choice: 1;
  #other choices, e.g., 32, 100
  division_const_tilde = 32 
)

###########################################
# Computation for poitwise sdf estimation
###########################################
#initialize
set.seed(cfg$seed)

# -----------------------
#Load data
# -----------------------
data = load_simulation_data(cfg, data_dir)

# true value of sdf at chosen frequency
sdf_true = sdf.true[length(sdf.true) * cfg$omega/pi]

#grid of privacy levels
alpha_grid =  exp(seq(log(0.01), log(100000), length.out = 100))

#theoretical truncation parameters
tau_NI_theo <- sqrt(56 * log(cfg$n)^(1 + cfg$delta))
tau_SI_theo <- sqrt(8 * log(cfg$n)^(1 + cfg$delta))

#list of truncation parameters to be used
#Note: labels in legend of the plots have to changed manually below,
#      in case other division factors, or more truncation parameters
#      are use than the default
tau_NI <- c(tau_NI_theo, tau_NI_theo / 6.5, tau_NI_theo / 20)
tau_SI <- c(tau_SI_theo, tau_SI_theo / 2.5, tau_SI_theo / 7)
tau_grid = cbind(tau_NI, tau_SI)

# -----------------------
# No privacy (no) baseline
# -----------------------
#number of covariance coef.s used for sdf estimation
K0 <- ceiling((1 / cfg$n)^(-1 / (2 * cfg$s + 1)))
mse_np <- numeric(cfg$N_MC)
for (i in seq_len(cfg$N_MC)) {
  X <- data[i, ]
  sdf_hat <- sdf.hat.global(X, K0, cfg$omega)
  mse_np[i] <- (sdf_hat - sdf_true)^2
}
MSE_np <- mean(mse_np)


# -----------------------
# Private estimators
# -----------------------
MSE_SI <- MSE_NI <- matrix(0, length(alpha_grid), length(tau_NI))

for (a in 1:length(alpha_grid)) {
  alpha = alpha_grid[a]
  
  for (t in 1:nrow(tau_grid)) {
    tau <- tau_grid[t, 1]
    tau_si <- tau_grid[t, 2]
    
    #number of covariance coef.s used for sdf estimation
    K_SI = ceiling(max(1 / cfg$n, tau_si^6 / (cfg$n * alpha^2))^(-1 / (2*cfg$s + 1)))
    K_NI = ceiling(max(1 / cfg$n, tau^4 / (cfg$n * alpha^4))^(-1 / (2*cfg$s + 1)))
    
    tau_tilde <- sqrt(1024 * tau_si^6 * (K_SI + 1)) / cfg$division_const_tilde
    
    mse_SI <- mse_NI <- numeric(cfg$N_MC)
    
    for (i in seq_len(cfg$N_MC)) {
      X <- data[i, ]
      
      # Sequential interactive privatization and estimation
      Z_SI <- priv.data.SI.sdf0(X, cfg$omega, alpha, tau_si, tau_tilde, K_SI)
      mse_SI[i] <- (sdf0.hat.SI(Z_SI) - sdf_true)^2
      
      # Non-interactive privatization and estimation
      Z_NI <- priv.data.NI(X, alpha, tau)
      sdf_NI <- sdf.hat.global.NI(Z_NI, K_NI, cfg$omega, alpha, tau)
      mse_NI[i] <- (sdf_NI - sdf_true)^2
      
    }
    MSE_SI[a, t] <- mean(mse_SI)
    MSE_NI[a, t] <- mean(mse_NI)
  }
}

results.SI <- results.NI <- list(
  target = "sdf",
  alpha_grid = alpha_grid,
  tau_grid = tau_grid,
  MSE_p = MSE_NI,
  MSE_np = MSE_np,
  sdf_true = sdf_true,
  cfg = cfg,
  typ = "NI"
)
results.SI$MSE_p =  MSE_SI
results.SI$typ = "SI"

#up to which alpha is has SI-based estimator smaller MSE?
#max(alpha_grid[results.SI$MSE_p[,1]<= results.NI$MSE_p[,1]])#theoretical tau
#max(alpha_grid[results.SI$MSE_p[,2]<= results.NI$MSE_p[,2]])#chosen tau

###############################################################
#Plots for pointwise sdf estimation: log(MSE) versus log(alpha)
###############################################################
#NI-plot
#manually adapt labels in case division factors for tau have been changed above
labels_NI = c(
  expression(tau == tau[theory]),
  bquote(tau == tau[theory] / 6.5),
  bquote(tau == tau[theory] / 20),
  "no privacy",
  expression(-4 * log(alpha))
)

#ylims = vector with lowest and highest value on y-axis,
#        in order to manually align axis of different plots
plot_mse_vs_alpha(
  results.NI,
  fig_dir,
  log_scale_y = TRUE,
  log_scale_x = TRUE,
  theo_rate = "NI",
  labels = labels_NI,
  ylims = c(-5,30)
)

#SI-plot
#manually adapt labels in case division factors for tau have been changed above
labels_SI = c(
  bquote(tau == tau[theory] * "," ~  ~ tilde(tau) == tilde(tau)[theory] / .(cfg$division_const_tilde)),
  bquote(tau == tau[theory] / 2.5 * "," ~  ~ tilde(tau) == tilde(tau)[theory] / .(cfg$division_const_tilde)),
  bquote(tau == tau[theory] / 7 * "," ~  ~ tilde(tau) == tilde(tau)[theory] / .(cfg$division_const_tilde)),
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

#save results
results_file =  
  sprintf(
    "%s%s_sdf0_n%d_omega%.2f_sd%.1f_s%d.RDATA",
    res_dir,
    cfg$experiment,
    cfg$n,
    cfg$omega,
    cfg$sd,
    cfg$s
  )

save(cfg, results.SI, results.NI, file = results_file)


# #######################################################################
# #Plot in paper
# #with two different division factors for tau.tilde
# #######################################################################
# # 1. run the above code with selected omega
# #    and cfg$division_const_tilde = 1 (example of first choice)
# #    optionally: to speed up set tau_NI <- c(tau_NI_theo) and tau_SI <- c(tau_SI_theo)
# #    above, but only in this run as other choices of tau not shown in plot below
# # 2. set
# results.SI_ti1= results.SI
# #3. run the above code with same omega
# #   and cfg$division_const_tilde = 32 (example of second choice)
# #   without running the line results.SI_ti1 = results.SI
# #4. execute
# results.SI_ti1$MSE_p=cbind(results.SI_ti1$MSE_p[,1],results.SI$MSE_p[,2:3])
# results.SI_ti1$tau_grid=rbind(results.SI_ti1$tau_grid[1,],results.SI$tau_grid[2:3,])
# # #5. Plot
# labels_SI = c(
#   bquote(tau == tau[theory] * "," ~~ tilde(tau) == tilde(tau)[theory]),
#   bquote(tau == tau[theory]/2.5 * "," ~~ tilde(tau) == tilde(tau)[theory]/ .(cfg$division_const_tilde)),
#   bquote(tau == tau[theory]/7 * "," ~~ tilde(tau) == tilde(tau)[theory]/ .(cfg$division_const_tilde)),
#   "no privacy",
#   expression(-2*log(alpha))
#   )
# #ylims = vector with lowest and highest value on y-axis,
# #        in order to manually align axis of different plots
# plot_mse_vs_alpha(
#   results.SI_ti1,
#   fig_dir,
#   log_scale_y = TRUE,
#   log_scale_x = TRUE,
#   theo_rate = "SI",
#   labels = labels_SI,
#   ylims = c(-5, 30)
# )
# 
# 
# results_file =  file.path(
#   "RESULTS",
#   sprintf(
#     "%s_sdf0_n%d_omega%.2f_sd%.1f_s%d.RDATA",
#     cfg$experiment, cfg$n, cfg$omega, cfg$sd, cfg$s
#   )
# )
# 
# save(cfg,results.SI_ti1,results.SI,results.NI,file=results_file)

source("data_simulation_functions.R")

######################
#CONFIGURATION
######################
#Before running: Create a folder DATA.

data_dir <- "../DATA/"

# Reproducibility
seed <- 42

# Model parameters
n <- 1000 #length of time series
sd <- 1.2 #standard deviation of Gaussian marginal distribution

# number of Monte Carlo iterations
N_MC <- 300


#Choose Example
ex.nb="EX1_"  #AR(0.8)-process
#ex.nb = "EX2_" #Hölder continuous sdf
#ex.nb="EX3_"  #acf with polynomial decay
#ex.nb="EX4_"  #white noise process


#Simulation of Data
set.seed(seed)

cat("Simulating data...\n")

data <- vector("list", N_MC)
if (ex.nb == "EX1_") {
  ex <- example1(n, N_MC, sd = sd)
  data <- ex$Y
  acf.true <- ex$acf
  sdf.true <- ex$sdf
}

if (ex.nb == "EX2_") {
  #gamma=1.7
  gamma = 0.8
  ex <- example2(n, N_MC, sd = sd, gamma = gamma)
  data <- ex$Y
  acf.true <- ex$acf
  sdf.true <- ex$sdf
}

if (ex.nb == "EX3_") {
  gamma = 5.1
  ex = example3(n, N_MC, sd = sd, gamma = gamma)
  data <- ex$Y
  acf.true <- ex$acf
  sdf.true <- ex$sdf
}

if (ex.nb == "EX4_") {
  acf.true = c(sd^2, rep(0, n - 1))
  sdf.true = acf2sdf(acf.true)
  data = matrix(mvrnorm(1, mu = numeric(n), Sigma = toeplitz(acf.true)), N_MC, n)
}


#file name
simdata_file <- file.path(sprintf("simdata_n%g_sd%.1f_MC%g_seed%g.RData", n, sd, N_MC, seed))
simdata_file <- paste0(data_dir, ex.nb, simdata_file)

#save data
save(data, acf.true, sdf.true, n, sd, N_MC, #gamma
     file = simdata_file)

cat("Saved to:\n", simdata_file, "\n")

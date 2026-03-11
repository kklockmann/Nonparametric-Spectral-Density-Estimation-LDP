#ylims = vector with lowest and highest value on y-axis,
#        in order to manually align axis of different plots
plot_mse_vs_alpha <- function(res,  
                                    fig_dir,
                                    log_scale_y = TRUE,
                                    log_scale_x = TRUE,
                                    theo_rate = NA,
                                    labels,
                                    ylims) {
  cfg <- res$cfg
  
  #Figure name
  if (res$target == "acf") {
    fname <- sprintf(
      "%s%s_acf_n%d_j%d_sd%.1f_%s.pdf",
      fig_dir,
      cfg$experiment,
      cfg$n,
      cfg$j,
      cfg$sd,
      res$typ
    )
  } else {
    fname <- sprintf(
      "%s%s_sdf0_n%d_omega%.2f_sd%.1f_s%d_%s.pdf",
      fig_dir,
      cfg$experiment,
      cfg$n,
      cfg$omega,
      cfg$sd,
      cfg$s,
      res$typ
    )
  }
  
  #set line typs
  nn <- dim(res$MSE_p)[2]
  if (nn < 4) {
    ltyps = c(seq(1, nn), 6, 5)
  } else{
    ltyps = c(setdiff(setdiff(seq(1, nn + 1), 5), 6), 6, 5)
  }
  #seq(1,t+1),
  
  alpha <- res$alpha_grid
  tau <- res$tau_grid[, 1] #first column is tau_grid
  
  for (t in 1:length(tau)) {
    MSE_np <- rep(res$MSE_np, length(alpha))
    MSE_p <- res$MSE_p[, t]
    
    if (log_scale_y) {
      MSE_p = log(MSE_p)
      MSE_np = log(MSE_np)
    }
    
    if (t == 1) {
      yvals <- c(MSE_p, MSE_np)
      ylim <- c(ylims[1], ylims[2]) #c(-5,30)
      xvals <- if (log_scale_x)
        log(alpha)
      else
        alpha
      
      pdf(fname, width = 10, height = 8)
      #par(mar = c(5, 6, 3, 2))
      par(mar = c(5, 6, 2.5, 2.5))
      
      plot(
        xvals,
        MSE_p,
        type = "l",
        lwd = 6,
        cex = 1.2,
        lty = t,
        ylim = ylim,
        ylab = "",
        xlab = "",
        cex.axis = 2.4
      )
      
    } else{
      lines(xvals, MSE_p, lwd = 4, lty = t)
    }
  }
  lines(xvals, MSE_np, lwd = 4, lty = 6)
  legend(
    "topright",
    legend = labels,
    lty = ltyps,
    lwd = 3,
    bty = "n",
    cex = 1.6,
    col = c(rep("black", nn + 1), "grey50")
  )
  
  mtext(
    if (log_scale_x)
      expression(Log(alpha))
    else
      expression(alpha),
    side = 1,
    line = 3.5,
    cex = 2.7
  )
  mtext(
    if (log_scale_y)
      "Estimated Log(MSE)"
    else
      "Estimated MSE",
    side = 2,
    line = 3.5,
    cex = 3.0
  )
  
  if (theo_rate == "NI" || theo_rate == "SI") {
    delta = cfg$delta
    if (theo_rate == "NI") {
      rate = log(log(cfg$n)^(2 + 2 * delta) / (cfg$n * alpha^4))
    } else{
      rate = log(log(cfg$n)^(4 + 4 * delta) / (cfg$n * alpha^2))
    }
    const = MSE_p[1] - rate[1]
    lines(xvals,
          rate,
          lwd = 3,
          lty = 5,
          col = "grey50")
  }
  
  dev.off()
  cat("Saved:", fname, "\n")
  
}

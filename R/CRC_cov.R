rm(list = ls())

packages <- c(
  "tidyverse",
  "dfoptim",
  "car",
  "plm",
  "gridExtra",
  "np",
  "AER",
  "Matrix",
  "nleqslv",
  "beepr"
)

pacman::p_load(packages, character.only = TRUE)

git_dir <- "~/Documents/GitHub/revisit_crc/"

data_simulate <- function(
  param, 
  seed, 
  rule = 1, # 1: take into account abs. advantage
            # 3. ignore abs. advantage
  N = 10000 # number of individual farmers
  ) {

  p           <- param[1]
  B_scale     <- param[2]
  B           <- param[3]
  sigma_N2    <- param[4]
  sigma_T2    <- param[5]
  sigma_NT    <- param[6]
  xi_scale    <- param[7]
  C_scale     <- param[8]
  C_level     <- param[9]
  Z_scale     <- param[10]
  rho         <- param[11]

  set.seed(seed)
  B_T <- B_scale + runif(1) / 10
  B_N <- B_T + B
  xi_N   <- matrix(rnorm(N * 2) * xi_scale, nrow = N)
  xi_T   <- matrix(rnorm(N * 2) * xi_scale, nrow = N)

  theta_sigma <- matrix(c(sigma_N2, sigma_NT, sigma_NT, sigma_T2), nrow = 2)
  lower_chol  <- t(chol(theta_sigma))
  theta_iNT   <- lower_chol %*% t(matrix(rnorm(N * 2), nrow = N))
  theta_iN    <- theta_iNT[1,]
  theta_iT    <- theta_iNT[2,]
  
  Z           <- matrix(rnorm(N * 2) * Z_scale, nrow = N)

  phi     <- (sigma_N2 - sigma_NT) / (sigma_NT - sigma_T2) - 1
  theta_i <- (sigma_NT - sigma_T2) / (sigma_N2 + sigma_T2 - 2 * sigma_NT) * (theta_iN - theta_iT)
  tau_i   <- theta_iT - theta_i

  C_N <- matrix(exp(cbind(runif(N), runif(N) + C_level) * C_scale), nrow = N)
  C_T <- matrix(exp(rnorm(2 * N) * C_scale), nrow = N)

  B_N_E_xi <- B_N + (xi_scale^2 / 2)
  B_T_E_xi <- B_T + (xi_scale^2 / 2)

  if (rule == 1) {
    Y_N <- exp(B_N_E_xi) * exp((phi + 1) * theta_i + rho * Z + tau_i)
    Y_T <- exp(B_T_E_xi) * exp(theta_i + rho * Z + tau_i)
    profit_N <- p * Y_N - C_N
    profit_T <- p * Y_T - C_T
    Y_N_post <- exp(B_N) * exp((phi + 1) * theta_i + tau_i + rho * Z + xi_N)
    Y_T_post <- exp(B_T) * exp(theta_i + tau_i + rho * Z + xi_T)
    h <- (profit_N >= profit_T)
  } else if (rule == 3) {
    Y_N <- exp(B_N_E_xi) * exp((phi + 1) * theta_i + rho * Z + tau_i)
    Y_T <- exp(B_T_E_xi) * exp(theta_i + rho * Z + tau_i)
    profit_N <- p * Y_N - C_N
    profit_T <- p * Y_T - C_T
    Y_N_post <- exp(B_N) * exp((phi + 1) * theta_i + tau_i + rho * Z + xi_N)
    Y_T_post <- exp(B_T) * exp(theta_i + tau_i + rho * Z + xi_T)
    h <-  (
      (p * exp(B_N_E_xi) * exp((phi + 1) * theta_i + rho * Z) - C_N) >= 
      (p * exp(B_T_E_xi) * exp(theta_i + rho * Z) - C_T)
      )
  }

  y <- h * log(Y_N_post) + (1 - h) * log(Y_T_post)
  profit <- h * profit_N + (1 - h) * profit_T
  profit_op <- (profit_N >= profit_T) * profit_N + (profit_N < profit_T) * profit_T

  return(list(y, profit, h, theta_i, profit_op, tau_i, Z))
}

min_dist_func <- function(
  str_param, 
  reduced_param, 
  h
  ) {
  
  lambda1 <- str_param[1]
  lambda2 <- str_param[2]
  lambda3 <- str_param[3]
  B       <- str_param[4]
  phi     <- str_param[5]

  h1 <- h[,1]
  h2 <- h[,2]
  h12 <- h1 * h2

  lambda0 <- - lambda1 * mean(h1) - lambda2 * mean(h2) - lambda3 * mean(h12)
  
  gamma1 <- reduced_param[1]
  gamma2 <- reduced_param[2]
  gamma3 <- reduced_param[3]
  gamma4 <- reduced_param[4]
  gamma5 <- reduced_param[5]
  gamma6 <- reduced_param[6]
  
  dist1 <- gamma1 - ((1 + phi) * lambda1 + B + phi * lambda0)
  dist2 <- gamma2 - lambda2
  dist3 <- gamma3 - ((1 + phi) * lambda3 + phi * lambda2)
  dist4 <- gamma4 - lambda1
  dist5 <- gamma5 - ((1 + phi) * lambda2 + B + phi * lambda0)
  dist6 <- gamma6 - ((1 + phi) * lambda3 + phi * lambda1)

  dist <- c(dist1, dist2, dist3, dist4, dist5, dist6)
  
  objective <- sqrt(dist %*% dist)
  return(objective)
  
}

structural_estimation <- function(
  y, 
  Z,
  h, 
  N = 10000,
  Z_ignore = FALSE
  ) {
  
  h1 <- h[,1]
  h2 <- h[,2]
  h12 <- h[,1] * h[,2]

  data <- data.frame(
    y = c(y[,1], y[,2]),
    Z = c(Z[,1], Z[,2]),
    h1 = rep(h1, 2), h2 = rep(h2, 2), h12 = rep(h12, 2), 
    time = c(rep(1, N), rep(2, N))
    )
  if (Z_ignore == FALSE) {
    res <- lm(
      y ~ (time == 2) + h1 + h2 + h12 + h1 * (time == 2) + h2 * (time == 2) + h12 * (time == 2) + Z, 
      data = data
      )
    test <- linearHypothesis(res, "h2TRUE = h1TRUE + time == 2TRUE:h1TRUE")$`Pr(>F)`[2]
    
    gamma1 <- res$coefficients[3]
    gamma2 <- res$coefficients[4]
    gamma3 <- res$coefficients[5]
    gamma4 <- res$coefficients[3] + res$coefficients[7]
    gamma5 <- res$coefficients[4] + res$coefficients[8]
    gamma6 <- res$coefficients[5] + res$coefficients[9]
    rho    <- res$coefficients[6]
  
    reduced_param <- c(gamma1, gamma2, gamma3, gamma4, gamma5, gamma6)
    
    str_param_init <- c(0.1, 0.1, 0.2, 0.7, -0.6)
  
    str_param <- nmk(
      par = str_param_init, 
      fn = min_dist_func, 
      reduced_param = reduced_param,
      h = h
      )$par
  
    str_param <- c(
      - str_param[1] * mean(h1) - str_param[2] * mean(h2) - str_param[3] * mean(h12), 
      str_param
      )
    return(list(str_param, test, rho))
  } else {
    res <- lm(
      y ~ (time == 2) + h1 + h2 + h12 + h1 * (time == 2) + h2 * (time == 2) + h12 * (time == 2),
      data = data
      )
    test <- linearHypothesis(res, "h2TRUE = h1TRUE + time == 2TRUE:h1TRUE")$`Pr(>F)`[2]
    
    gamma1 <- res$coefficients[3]
    gamma2 <- res$coefficients[4]
    gamma3 <- res$coefficients[5]
    gamma4 <- res$coefficients[3] + res$coefficients[6]
    gamma5 <- res$coefficients[4] + res$coefficients[7]
    gamma6 <- res$coefficients[5] + res$coefficients[8]
  
    reduced_param <- c(gamma1, gamma2, gamma3, gamma4, gamma5, gamma6)
    
    str_param_init <- c(0.1, 0.1, 0.2, 0.7, -0.6)
  
    str_param <- nmk(
      par = str_param_init, 
      fn = min_dist_func, 
      reduced_param = reduced_param,
      h = h
      )$par
  
    str_param <- c(
      - str_param[1] * mean(h1) - str_param[2] * mean(h2) - str_param[3] * mean(h12), 
      str_param
      )
    return(list(str_param, test))
    
  }
}

monte_carlo_simulation <- function(
  param, 
  numSim, 
  rule = 1,
  Z_ignore = FALSE
  ) {
  
  monte_carlo_output <- c()
  for (b in 1:numSim) {
    seed <- b
    data <- data_simulate(param, seed, rule)

    p           <- param[1]
    B_scale     <- param[2]
    B           <- param[3]
    sigma_N2    <- param[4]
    sigma_T2    <- param[5]
    sigma_NT    <- param[6]
    xi_scale    <- param[7]
    C_scale     <- param[8]
    C_level     <- param[9]
    Z_scale     <- param[10]
    rho         <- param[11]
    phi         <- (sigma_N2 - sigma_NT) / (sigma_NT - sigma_T2) - 1

    y         <- data[[1]]
    profit    <- data[[2]]
    h         <- data[[3]]
    theta_i   <- data[[4]]
    profit_op <- data[[5]]
    Z         <- data[[7]]
  
    estimation <- structural_estimation(y, Z, h, Z_ignore = Z_ignore)
    str_param <- estimation[[1]]
    p_value   <- estimation[[2]]
    if (Z_ignore == FALSE) {
      rho_est   <- estimation[[3]]
    }
    lambda0_est  <- str_param[1]
    lambda1_est  <- str_param[2]
    lambda2_est  <- str_param[3]
    lambda12_est <- str_param[4]
    B_est        <- str_param[5]
    phi_est      <- str_param[6]
  
    theta_i0_true  <- mean(theta_i[h[,1] == 0 & h[,2] == 0])
    theta_i1_true  <- mean(theta_i[h[,1] == 1 & h[,2] == 0])
    theta_i2_true  <- mean(theta_i[h[,1] == 0 & h[,2] == 1])
    theta_i12_true <- mean(theta_i[h[,1] == 1 & h[,2] == 1])
  
    theta_i0_model  <- lambda0_est
    theta_i1_model  <- lambda0_est + lambda1_est
    theta_i2_model  <- lambda0_est + lambda2_est
    theta_i12_model <- lambda0_est + lambda1_est + lambda2_est + lambda12_est
  
    theta_i0_diff  <- (theta_i0_model  - theta_i0_true) 
    theta_i1_diff  <- (theta_i1_model  - theta_i1_true) 
    theta_i2_diff  <- (theta_i2_model  - theta_i2_true) 
    theta_i12_diff <- (theta_i12_model - theta_i12_true)

    return_i0_diff  <- (
      (B_est + phi_est * theta_i0_model)  - 
        (B + phi * theta_i0_true)
      ) 
    return_i1_diff  <- (
      (B_est + phi_est * theta_i1_model) - 
        (B + phi * theta_i1_true)
      ) 
    return_i2_diff  <- (
      (B_est + phi_est * theta_i2_model) -
        (B + phi * theta_i2_true)
      ) 
    return_i12_diff <- (
      (B_est + phi_est * theta_i12_model) - 
        (B + phi * theta_i12_true)
      )
  
    percent_loss <- mean((profit - profit_op) / abs(profit_op))
  
    output <- c(
      theta_i0_diff, theta_i1_diff, theta_i2_diff, theta_i12_diff, 
      B_est, phi_est, p_value, percent_loss,
      return_i0_diff, return_i1_diff, return_i2_diff, return_i12_diff,
      theta_i0_model, theta_i1_model, theta_i2_model, theta_i12_model
      )
    if (Z_ignore == FALSE) {
      output <- c(output, rho_est)
    }
  
    monte_carlo_output <- rbind(monte_carlo_output, output)
  }

  return(monte_carlo_output)
}

make_fig <- function(est_return, label, x_lim, binwidth = 0.025) {
  print(
    ggplot() + 
      geom_histogram(aes(x = est_return), binwidth = binwidth) +
      coord_cartesian(xlim = x_lim) +
      xlab(label) +
      ylab("Count") +
      geom_vline(xintercept = 0, color = "red", linewidth = 1.0, alpha = 0.5) +
      theme_minimal()
  )
}

# Monte Carlo Simulation
numSim <- 1000

p           <- 1
B_scale     <- 0.3
B           <- 1.0
sigma_N2    <- 3.5
sigma_T2    <- 10.5
sigma_NT    <- 6.0
xi_scale    <- 0.1
C_scale     <- 1.1
C_level     <- 1.6
Z_scale     <- 1.0
rho         <- 2.0
phi         <- (sigma_N2 - sigma_NT) / (sigma_NT - sigma_T2) - 1

param <- c(
  p, B_scale, B,
  sigma_N2, sigma_T2, sigma_NT, xi_scale,
  C_scale, C_level, Z_scale, rho
  )

# Monte Carlo Simulations =====================
# When maximizing expected profits
monte_carlo_covar_pm <- monte_carlo_simulation(
  param, numSim, rule = 1
  )
monte_carlo_covar_pm[,7] %>% plot

# Figures when maximizing expected profits ===================
# return ------------
x_lim <- c(-0.5, 0.5)
covar_pm_return_0 <- make_fig(monte_carlo_covar_pm[,9], "Never adopters", x_lim)
covar_pm_return_1 <- make_fig(monte_carlo_covar_pm[,10], "Early adopters", x_lim)
covar_pm_return_2 <- make_fig(monte_carlo_covar_pm[,11], "Late adopters", x_lim)
covar_pm_return_12 <- make_fig(monte_carlo_covar_pm[,12], "Always adopters", x_lim)
ggsave(
  filename = file.path(
    git_dir,
    "Figures/covar_pm_return.pdf"
  ),
  plot = grid.arrange(
    covar_pm_return_0, covar_pm_return_1, covar_pm_return_2, covar_pm_return_12, 
    ncol = 2
    ),
  height = 7,
  width = 7
  )
 
# Figures when not maximizing expected profits ===================
# When not maximizing expected profits
monte_carlo_covar_nopm <- monte_carlo_simulation(
  param, numSim, rule = 3
  )
monte_carlo_covar_nopm[,17] %>% plot
monte_carlo_covar_nopm[,7] %>% plot
(monte_carlo_covar_nopm[,7] < 0.1) %>% mean

# return ------------
x_lim <- c(-2.0, 2.0)
covar_nopm_return_0 <- make_fig(monte_carlo_covar_nopm[,9], "Never adopters", x_lim, binwidth = 0.2)
covar_nopm_return_1 <- make_fig(monte_carlo_covar_nopm[,10], "Early adopters", x_lim, binwidth = 0.2)
covar_nopm_return_2 <- make_fig(monte_carlo_covar_nopm[,11], "Late adopters", x_lim, binwidth = 0.2)
covar_nopm_return_12 <- make_fig(monte_carlo_covar_nopm[,12], "Always adopters", x_lim, binwidth = 0.2)
ggsave(
  filename = file.path(
    git_dir,
    "Figures/covar_nopm_return.pdf"
  ),
  plot = grid.arrange(
    covar_nopm_return_0,
    covar_nopm_return_1,
    covar_nopm_return_2,
    covar_nopm_return_12,
    ncol=2
    ),
  height = 7,
  width = 7
  )
 
# rho
covar_nopm_rho <- ggplot() +
  geom_histogram(aes(x=monte_carlo_covar_nopm[,17]), binwidth=0.01) + 
  coord_cartesian(xlim = c(1.9, 2.4)) +
  xlab("Estimated rho") +
  ylab("Count") +
  geom_vline(xintercept = rho, color = "red", linewidth = 1.0, alpha = 0.5) +
  theme_minimal()
ggsave(
  filename = file.path(
    git_dir,
    "Figures/covar_nopm_rho.pdf"
  ),
  plot = covar_nopm_rho,
  height = 7,
  width = 7
  )



# Figures when not maximizing expected profits ===================
# When not maximizing expected profits
monte_carlo_covar_nopm_Z_ignore <- monte_carlo_simulation(
  param, numSim, rule = 3, Z_ignore = TRUE
  )
monte_carlo_covar_nopm_Z_ignore[,7] %>% plot
(monte_carlo_covar_nopm_Z_ignore[,7] < 0.1) %>% mean

# return ------------
x_lim <- c(-3.5, 3.5)
covar_nopm_return_0_Z_ignore <- make_fig(
  monte_carlo_covar_nopm_Z_ignore[monte_carlo_covar_nopm_Z_ignore[,7] < 0.1,9],
  "Never adopters", x_lim, binwidth = 0.2 
  )
covar_nopm_return_1_Z_ignore <- make_fig(
  monte_carlo_covar_nopm_Z_ignore[monte_carlo_covar_nopm_Z_ignore[,7] < 0.1,10],
  "Early adopters", x_lim, binwidth = 0.2
  )
covar_nopm_return_2_Z_ignore <- make_fig(
  monte_carlo_covar_nopm_Z_ignore[monte_carlo_covar_nopm_Z_ignore[,7] < 0.1, 11],
  "Late adopters", x_lim, binwidth = 0.2
  )
covar_nopm_return_12_Z_ignore <- make_fig(
  monte_carlo_covar_nopm_Z_ignore[monte_carlo_covar_nopm_Z_ignore[,7] < 0.1, 12], 
  "Always adopters", x_lim, binwidth = 0.2
  )
ggsave(
  filename = file.path(
    git_dir,
    "Figures/covar_nopm_return_Z_ignore.pdf"
  ),
  plot = grid.arrange(
    covar_nopm_return_0_Z_ignore,
    covar_nopm_return_1_Z_ignore,
    covar_nopm_return_2_Z_ignore,
    covar_nopm_return_12_Z_ignore,
    ncol=2
    ),
  height = 7,
  width = 7
  )
 


x <- c(22, 27)               
n_total <- c(319, 319)      
y_infected <- round(c(20.61, 17.55))

x0 <- 14                              
n_power <- 2                         
f_model <- function(x, A, B, C) {
  d <- x - x0
  B + (A * d^n_power - B) * exp(-C * d)
}

neg_log_likelihood <- function(par) {
  A <- par[1]
  B <- par[2]
  C <- par[3]
  
  p <- f_model(x, A, B, C)
  p <- pmax(pmin(p, 1 - 1e-6), 1e-6)  
  
  -sum(y_infected * log(p) + (n_total - y_infected) * log(1 - p))
}

result <- optim(
  par = c(0.01, 0.1, 0.1),
  fn = neg_log_likelihood,
  method = "BFGS",
  hessian = TRUE
)

par_hat <- result$par
cov_matrix <- solve(result$hessian)

ages <- 20:27
n_new <- 319  

pred_df <- data.frame(
  Age = ages,
  Predicted_Prevalence = NA,
  CI_Lower = NA,
  CI_Upper = NA,
  Predicted_Infections = NA
)

for (i in seq_along(ages)) {
  x_new <- ages[i]
  d <- x_new - x0
  A <- par_hat[1]; B <- par_hat[2]; C <- par_hat[3]
  
  p_hat <- f_model(x_new, A, B, C)
  p_hat <- pmax(pmin(p_hat, 1 - 1e-6), 1e-6)
  
  grad <- numeric(3)
  grad[1] <- d^n_power * exp(-C * d)
  grad[2] <- 1 - exp(-C * d)
  grad[3] <- -(A * d^n_power - B) * d * exp(-C * d)
  
  var_pred <- t(grad) %*% cov_matrix %*% grad
  se_pred <- sqrt(var_pred)
  
  # CI
  ci_lower <- p_hat - 1.96 * se_pred
  ci_upper <- p_hat + 1.96 * se_pred
  ci_lower <- pmax(ci_lower, 0)
  ci_upper <- pmin(ci_upper, 1)
  
  pred_df$Predicted_Prevalence[i] <- round(p_hat, 4)
  pred_df$CI_Lower[i] <- round(ci_lower, 4)
  pred_df$CI_Upper[i] <- round(ci_upper, 4)
  pred_df$Predicted_Infections[i] <- round(p_hat * n_new, 2)
}

print(pred_df)


age_seq <- 20:29
newdata <- data.frame(x = age_seq, n = 1)
posterior_preds <- posterior_epred(fit_brm, newdata = newdata)
f_mean <- apply(posterior_preds, 2, mean)
f_lower <- apply(posterior_preds, 2, quantile, probs = 0.025)
f_upper <- apply(posterior_preds, 2, quantile, probs = 0.975)
result_df <- data.frame(
  Age = age_seq,
  Predicted_Prevalence = round(f_mean, 4),
  Lower_95_CI = round(f_lower, 4),
  Upper_95_CI = round(f_upper, 4)
)
print(result_df, row.names = FALSE)

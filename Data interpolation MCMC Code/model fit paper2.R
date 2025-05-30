x <- c(16.5, 22, 27, 34.5)               
n_total <- c(740, 445, 414, 903)      
y_infected <- c(26.64, 72.09, 42.64, 51.47)  
y_infected <- round(c(26.64, 72.09, 42.64, 51.47))

x0 <- 12                             
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
  method = "BFGS"
)


curve_x <- seq(14, 35, length.out = 200)
fitted_p <- f_model(curve_x, result$par[1], result$par[2], result$par[3])

plot(x, y_infected / n_total, col = "red", pch = 19,
     ylim = c(0, max(y_infected / n_total, fitted_p) + 0.05),
     xlab = "Age", ylab = "HPV prevalence",
     main = "Fitted curve")
lines(curve_x, fitted_p, col = "blue", lwd = 2)
legend("topright", legend = c("Observed", "Fitted"),
       col = c("red", "blue"), pch = c(19, NA), lty = c(NA, 1))








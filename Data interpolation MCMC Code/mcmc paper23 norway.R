library(brms)
df <- data.frame(
  x = c(22, 38.5),
  y = c(228, 75),        
  n = c(1011, 1217)     
)

df$y <- round(df$y)
formula_nl <- bf(
  y | trials(n) ~ B + (A * (x - 14)^2 - B) * exp(-C * (x - 14)),
  A + B + C ~ 1,
  nl = TRUE
)
priors <- c(
  prior(normal(0.0757, 0.01), nlpar = "A"),
  prior(normal(-0.1027, 0.01), nlpar = "B"),
  prior(normal(0.0963, 0.01), nlpar = "C")
)

fit_brm <- brm(
  formula = formula_nl,
  data = df,
  family = binomial(link = "identity"),  
  prior = priors,
  chains = 4,
  iter = 4000,
  seed = 123,
  warmup = 1000
)

summary(fit_brm)

newdata <- data.frame(x = 16, n = 1)
pred <- posterior_epred(fit_brm, newdata = newdata)
ci <- quantile(pred, probs = c(0.025, 0.975))

plot(fit_brm)
mcmc_plot(fit_brm, type = "dens_overlay")
bayesplot::mcmc_acf()

age_seq <- 18:40
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



plot(age_seq, f_mean, type = "l", lwd = 2, col = "black",
     ylim = c(0, max(f_upper) + 0.05),
     xlab = "Age", ylab = "Predicted HPV Prevalence",
     main = "HPV Prevalence Prediction")


polygon(c(age_seq, rev(age_seq)),
        c(f_lower, rev(f_upper)),
        col = rgb(0.7, 0.7, 0.7, 0.4), border = NA)
lines(age_seq, f_mean, col = "black", lwd = 2)

points(df$x, df$y / df$n, col = "blue", pch = 19)
legend("topright", legend = c("Observed", "Predicted", "95% CI"),
       col = c("black", "blue", rgb(0.7, 0.7, 0.7, 0.4)),
       pch = c(19, NA, NA), lty = c(NA, 1, NA), lwd = c(NA, 2, NA),
       pt.cex = 1.2, bty = "n")





library(brms)
df <- data.frame(
  x = c(16.5, 22, 27, 34.5),
  y = c(26.64, 72.09, 42.64, 51.47) ,
  n = c(740, 445, 414, 903)
)
df$y <- round(df$y)
formula_nl <- bf(
  y | trials(n) ~ B + (A * (x - 14)^2 - B) * exp(-C * (x - 14)),
  A + B + C ~ 1,
  nl = TRUE
)
priors <- c(
  prior(normal(0.01, 0.05), nlpar = "A"),        
  prior(normal(0.05, 0.05), nlpar = "B"),        
  prior(normal(0.2, 0.1),  nlpar = "C")          
)

fit_brm <- brm(
  formula = formula_nl,
  data = df,
  family = binomial(link = "identity"),  
  prior = priors,
  chains = 4,
  iter = 4000,
  seed = 123,
  warmup = 1000,
  control = list(adapt_delta = 0.99)
  
)

summary(fit_brm)

newdata <- data.frame(x = 16, n = 1)
pred <- posterior_epred(fit_brm, newdata = newdata)
ci <- quantile(pred, probs = c(0.025, 0.975))

plot(fit_brm)
mcmc_plot(fit_brm, type = "dens_overlay")
bayesplot::mcmc_acf()

 

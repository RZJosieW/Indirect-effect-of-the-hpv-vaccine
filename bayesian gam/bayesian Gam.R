install.packages("brms")
install.packages("devtools")
library("brms")
library("mgcv")
library("ggplot2")
library("schoenberg")
library("rstan")
devtools::install_github('gavinsimpson/schoenberg')
library(readxl)
hpvdata <- read_excel("/Users/V/Desktop/hpvfinal.xlsx")
head(hpvdata)
colnames(hpvdata)[colnames(hpvdata) == "Pre HPV-16, -18"] <- "pre"
colnames(hpvdata)[colnames(hpvdata) == "Post HPV-16, -18"] <- "post"
colnames(hpvdata)[colnames(hpvdata) == "vaccine coverage"] <- "coverage"
library(ggplot2)
library(gratia)  
library(mgcv)
hpv_clean <- hpvdata[, c("indirect", "coverage", "agemid", "year", "id", "type","sdpost")]
hpv_clean <- hpv_clean[-213, ]






prior.list <- get_prior(
  indirect ~ 
    s(coverage, k = 5) + 
    s(agemid, k = 5) + 
    s(year, k = 5),
  data = hpv_clean,
  family = gaussian()
)

ptm <- proc.time()

model.bayes.gam <- brm(bf(indirect ~  s(coverage, k = 5) + s(agemid, k = 5) + s(year, k = 5)),
                       data = hpv_clean, 
                       family = gaussian(),
                       prior = prior.list,
                    
                       iter = 4000, 
                       warmup = 1000, 
                      
                       refresh = 0,
                       control = list(adapt_delta = 0.99))

proc.time() - ptm
summary(model.bayes.gam)

plot(conditional_smooths(model.bayes.gam), ask = FALSE)
library(bayesplot)
posterior <- as_draws_df(model.bayes.gam)
posterior <- as_draws_df(model.bayes.gam)
library(bayesplot)
mcmc_areas(posterior, pars = c("b_Intercept", "sds_scoverage_1", "sds_sagemid_1", "sds_syear_1"))
pp_check(model.bayes.gam)








# t distrubution
prior.list <- get_prior(
  indirect ~ 
    s(coverage, k = 5) + 
    s(agemid, k = 5) + 
    s(year, k = 5),
  data = hpv_clean,
  family = gaussian()
)


priors <- c(
  prior(normal(0, 10), class = "Intercept"),
  prior(exponential(2), class = "sds"),    
  prior(cauchy(0, 2.5), class = "sigma")
)

ptm <- proc.time()
model.bayes.gam1 <- brm(bf(indirect ~  s(coverage, k = 5) + s(agemid, k = 5) + s(year, k = 5)),
                       data = hpv_clean, 
                       family = student(),
                       prior = priors,
                       
                       iter = 4000, 
                       warmup = 1000, 
                       
                       refresh = 0,
                       control = list(adapt_delta = 0.99))
proc.time() - ptm
summary(model.bayes.gam1)

plot(conditional_smooths(model.bayes.gam), ask = FALSE)
library(bayesplot)
posterior <- as_draws_df(model.bayes.gam)
posterior <- as_draws_df(model.bayes.gam)
library(bayesplot)
mcmc_areas(posterior, pars = c("b_Intercept", "sds_scoverage_1", "sds_sagemid_1", "sds_syear_1"))
pp_check(model.bayes.gam)






smooth_plot <- conditional_smooths(model.bayes.gam)
smooth_plot$`mu: s(coverage)`




# change k and exp distribution

priors <- c(
  prior(normal(0, 10), class = "Intercept"),
  prior(exponential(1), class = "sds"),    
  prior(cauchy(0, 2.5), class = "sigma")
)

model.bayes.gam2 <- brm(
  bf(indirect ~ 
       s(coverage, k = 8) +  
       s(agemid, k = 8) + 
       s(year, k = 5)), 
  data = hpv_clean,
  family = student(),
  prior = priors,
  chains = 4, iter = 4000, warmup = 1000,
  control = list(adapt_delta = 0.99)
)

summary(model.bayes.gam2)
library(bayesplot)
posterior2 <- as_draws_df(model.bayes.gam2)
mcmc_areas(posterior2, pars = c("b_Intercept", "sds_scoverage_1", "sds_sagemid_1", "sds_syear_1"))
pp_check(model.bayes.gam2)

smooth_plot2 <- conditional_smooths(model.bayes.gam2)
smooth_plot2$`s(coverage)`
smooth_plot2$`s(agemid)`
smooth_plot2$`s(year)`

smooth_plot2 <- conditional_smooths(model.bayes.gam2)
names(smooth_plot2)
smooth_plot2$`mu: s(coverage,k=8)`
smooth_plot2$`mu: s(agemid,k=8)`
smooth_plot2$`mu: s(year,k=5)`




coverage_smooth <- as.data.frame(smooth_plot2$`mu: s(coverage,k=8)`)
agemid_smooth <- as.data.frame(smooth_plot2$`mu: s(agemid,k=8)`)
year_smooth <- as.data.frame(smooth_plot2$`mu: s(year,k=5)`)

library(ggplot2)

ggplot(coverage_smooth, aes(x = coverage, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "steelblue3") +   
  geom_line(color = "navyblue", size = 1.2) +   
  labs(
    title = "Posterior Estimates of the Indirect Effect by Vaccine Coverage",
    x = "Vaccine Coverage",
    y = "Indirect Effect"
  ) +
  theme_minimal(base_size = 14)



ggplot(agemid_smooth, aes(x = agemid, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "seagreen3") +
  geom_line(color = "seagreen4", size = 1.2) +
  labs(
    title = "Posterior Estimates of the Indirect Effect by Age",
    x = "Age",
    y = "Indirect Effect"
  ) +
  theme_minimal(base_size = 14)



ggplot(year_smooth, aes(x = year, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "lightcoral") +
  geom_line(color = "indianred3", size = 1.2) +
  labs(
    title = "Posterior Estimates of the Indirect Effect by Years Since Vaccine Introduction",
    x = "Years Since Vaccine Introduction",
    y = "Indirect Effect"
  ) +
  theme_minimal(base_size = 14)




library(brms)

bayes_R2(model.bayes.gam2)
fitted_values <- fitted(model.bayes.gam2)[, "Estimate"]
actual_values <- hpv_clean$indirect
mse <- mean((fitted_values - actual_values)^2)
mae <- mean(abs(fitted_values - actual_values))
mse
mae


library(brms)
 post_pred <- posterior_predict(model.bayes.gam2)  
 y_true <- hpv_clean$indirect  
 n_iter <- nrow(post_pred)     
 
 mse_samples <- rep(NA, n_iter)
 mae_samples <- rep(NA, n_iter)
 
 for (i in 1:n_iter) {
       y_pred <- post_pred[i, ]    
       mse_samples[i] <- mean((y_true - y_pred)^2)  
       mae_samples[i] <- mean(abs(y_true - y_pred))  
     }
mse_point <- mean(mse_samples)
mae_point <- mean(mae_samples)
 

 quantile(mse_samples, probs = c(0.025, 0.975))  

 quantile(mae_samples, probs = c(0.025, 0.975))  
post_pred <- posterior_predict(model.bayes.gam2)  










library(brms)
library(caret)



set.seed(123)

folds <- createFolds(hpv_clean$indirect, k = 5)

mse_draws_all <- c()
mae_draws_all <- c()

for (i in 1:5) {
  test_idx <- folds[[i]]
  train_data <- hpv_clean[-test_idx, ]
  test_data  <- hpv_clean[test_idx, ]
  y_test <- test_data$indirect
  
  fit_cv <- brm(
    bf(indirect ~ s(coverage, k = 8) + s(agemid, k = 8) + s(year, k = 5)),
    data = train_data,
    family = student(),
    prior = priors,
    chains = 2, iter = 2000, warmup = 1000, refresh = 0,
    control = list(adapt_delta = 0.99)
  )
  bayes_R2(fit_cv)
  
  post_preds <- posterior_predict(fit_cv, newdata = test_data)  # [draws × test set size]
  
  for (j in 1:nrow(post_preds)) {
    y_pred <- post_preds[j, ]
    
    # MSE / MAE
    mse_draws_all <- c(mse_draws_all, mean((y_test - y_pred)^2))
    mae_draws_all <- c(mae_draws_all, mean(abs(y_test - y_pred)))
    
    # R²
    ss_res <- sum((y_test - y_pred)^2)
    ss_tot <- sum((y_test - mean(y_test))^2)
    r2 <- 1 - ss_res / ss_tot
    r2_draws_all <- c(r2_draws_all, r2)
  }
}


mse_cv_point <- mean(mse_draws_all)
mae_cv_point <- mean(mae_draws_all)
r2_cv_point  <- mean(r2_draws_all)

mse_cv_ci <- quantile(mse_draws_all, probs = c(0.025, 0.975))
mae_cv_ci <- quantile(mae_draws_all, probs = c(0.025, 0.975))
r2_cv_ci  <- quantile(r2_draws_all,  probs = c(0.025, 0.975))



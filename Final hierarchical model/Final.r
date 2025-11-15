library(readxl)
library(dplyr)
library(brms)
library(ggplot2)

hpvdata <- read_excel("/Users/V/Desktop/hpcdatafinalse.xlsx")
colnames(hpvdata)[colnames(hpvdata) == "vaccine coverage"] <- "coverage"


hpvdata$id <- factor(hpvdata$id)

hpvdata$countries <- factor(hpvdata$countries, levels = 1:11)
hpvdata$se <- as.numeric(hpvdata$SE)
hpvdata$type <- factor(as.integer(hpvdata$type),
                       levels = c(1, 2, 3),
                       labels = c("routine", "catchup", "unvaccinated"))
hpvdata$countries <- factor(hpvdata$countries,
                            levels = sort(unique(hpvdata$countries)))

# hpv list
hpv <- hpvdata %>%
  filter(!is.na(indirect), !is.na(coverage), !is.na(year),
         !is.na(agemid),  !is.na(type),     !is.na(countries),
         !is.na(se),
         !is.na(id)) %>%  
  mutate(
    coverage  = as.numeric(coverage),
    year      = as.numeric(year),
    agemid    = as.numeric(agemid),
    se        = as.numeric(se),
    type      = factor(type,      levels = c("routine","catchup","unvaccinated")),
    countries = factor(countries, levels = levels(hpvdata$countries)),
    id        = factor(id)  #
  ) %>%
  droplevels()

stopifnot(all(c("indirect","coverage","year","agemid","type","countries", "id") %in% names(hpv))) 



pri_B <- c(
  prior(student_t(3, 0, 5), class = "Intercept"),
  prior(student_t(3, 0, 5), class = "b"),
  prior(gamma(2, 0.1),       class = "nu"),    
  prior(exponential(3),        class = "sds"),   
  
  prior(exponential(1.5),    class = "sd", group = "countries:type"),
  prior(exponential(1.5),    class = "sd", group = "id")
)

# ---  
form_B <- bf(
  indirect | se(se, sigma = TRUE) ~ 1 + type +
    s(coverage, bs = "tp", k = 8) +
    s(year,     bs = "tp", k = 6) +
    s(agemid,   bs = "ts", k = 6) +
    s(coverage, by = type, bs = "fs", k = 12) +
    s(year,     by = type, bs = "fs", k = 8) +
    s(agemid,   by = type, bs = "ts", k = 8) +
    (1 | countries:type) +
    (1 | id)  
)






# --- fit ---
fit_B <- brm(
  form_B, data = hpv,
  family = student(),
  prior = pri_B,
  chains = 4, cores = 4, iter = 4000, warmup = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 13)
)



summary(fit_B)
print(brms::pp_check(fit_B, nsamples = 50))
brms::check_hmc_diagnostics(fit_B)
res <- residuals(fit_B, type = "pearson")
plot(fitted(fit_B), res)
abline(h=0, lty=2)

loo_pit_vals <- brms::loo_pit(fit_B)
print(bayesplot::ppc_loo_pit_qq(loo_pit_vals))
plot(fit_B)



# --- Diagnosis ---
summary(fit_B)
print(brms::pp_check(fit_B, nsamples = 100))
brms::check_hmc_diagnostics(fit_B)
res <- residuals(fit_B, type = "pearson")
plot(fitted(fit_B), res)
abline(h=0, lty=2)



library(loo)
loo_pit_vals <- loo_pit(fit_B) 
loo_pit_vals <- loo::loo_pit(fit_B)
print(bayesplot::ppc_loo_pit_qq(loo_pit_vals))


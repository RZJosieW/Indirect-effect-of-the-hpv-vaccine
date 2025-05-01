library(haven)
library(haven)
library(dplyr)
library(survey)


path <- "/Users/V/Desktop/uscoverage/"

# 
data_0910 <- inner_join(
  read_xpt(paste0(path, "0910.XPT")) %>% mutate(cycle = "2009-2010"),
  read_xpt(paste0(path, "D0910.XPT")) %>% mutate(cycle = "2009-2010"),
  by = c("SEQN", "cycle")
)

data_1112 <- inner_join(
  read_xpt(paste0(path, "1112.XPT")) %>% mutate(cycle = "2011-2012"),
  read_xpt(paste0(path, "D1112.XPT")) %>% mutate(cycle = "2011-2012"),
  by = c("SEQN", "cycle")
)

combined_0912 <- bind_rows(data_0910, data_1112)

filtered_0912 <- combined_0912 %>%
  filter(RIAGENDR == 2) %>%
  select(SEQN, cycle, RIDAGEYR, RIAGENDR, IMQ040, WTINT2YR, SDMVPSU, SDMVSTRA) %>%
  mutate(WTINT8YR = WTINT2YR / 2)  

filtered_0912_valid <- filtered_0912 %>%
  filter(
    RIDAGEYR >= 9 & RIDAGEYR <= 59,
    IMQ040 %in% c(1, 2)
  )

svy_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTINT8YR,
  nest = TRUE,
  data = filtered_0912_valid
)

weighted_coverage <- svyby(
  ~I(IMQ040 == 1),
  ~RIDAGEYR,
  svy_design,
  svymean,
  na.rm = TRUE
)

weighted_coverage_clean <- weighted_coverage %>%
  select(RIDAGEYR, `I(IMQ040 == 1)TRUE`, `se.I(IMQ040 == 1)TRUE`) %>%
  rename(
    vaccine_coverage = `I(IMQ040 == 1)TRUE`,
    se_vaccine_coverage = `se.I(IMQ040 == 1)TRUE`
  ) %>%
  mutate(
    lower_CI = vaccine_coverage - 1.96 * se_vaccine_coverage,
    upper_CI = vaccine_coverage + 1.96 * se_vaccine_coverage
  )

print(weighted_coverage_clean)

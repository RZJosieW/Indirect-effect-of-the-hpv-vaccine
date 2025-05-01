library(haven)
library(haven)
library(dplyr)
library(survey)

path <- "/Users/V/Desktop/uscoverage/"

data_0910 <- inner_join(
  read_xpt(paste0(path, "0910.XPT")) %>% mutate(cycle = "2009-2010"),
  read_xpt(paste0(path, "D0910.XPT")) %>% mutate(cycle = "2009-2010"),
  by = c("SEQN", "cycle")
)

filtered_0910 <- data_0910 %>%
  filter(RIAGENDR == 2) %>%
  select(SEQN, cycle, RIDAGEYR, RIAGENDR, IMQ040, IMQ045,
         WTINT2YR, SDMVPSU, SDMVSTRA)
data_0708 <- N %>%
  mutate(cycle = "2007-2008") %>%
  filter(RIAGENDR == 2) %>%
  select(SEQN, cycle, RIDAGEYR, RIAGENDR, IMQ040, IMQ045,
         WTINT2YR, SDMVPSU, SDMVSTRA)

full_data0710 <- bind_rows(data_0708, filtered_0910)


full_data0710 <- full_data0710 %>%
  mutate(WTINT8YR = WTINT2YR / 2)  


filtered_0710_valid <- full_data0710 %>%
  filter(
    RIDAGEYR >= 9 & RIDAGEYR <= 59,
    IMQ040 %in% c(1, 2)  
  )

svy_design_0710 <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTINT8YR,
  nest = TRUE,
  data = filtered_0710_valid
)

weighted_coverage_0710 <- svyby(
  ~I(IMQ040 == 1),
  ~RIDAGEYR,
  svy_design_0710,
  svymean,
  na.rm = TRUE
)

weighted_coverage_clean_0710 <- weighted_coverage_0710 %>%
  select(RIDAGEYR, `I(IMQ040 == 1)TRUE`, `se.I(IMQ040 == 1)TRUE`) %>%
  rename(
    vaccine_coverage = `I(IMQ040 == 1)TRUE`,
    se_vaccine_coverage = `se.I(IMQ040 == 1)TRUE`
  ) %>%
  mutate(
    lower_CI = vaccine_coverage - 1.96 * se_vaccine_coverage,
    upper_CI = vaccine_coverage + 1.96 * se_vaccine_coverage
  )

print(weighted_coverage_clean_0710)

library(haven)
library(dplyr)
library(survey)

path <- "/Users//Desktop/uscoverage/"

data_1112 <- inner_join(
  read_xpt(paste0(path, "1112.XPT")) %>% mutate(cycle = "2011-2012"),
  read_xpt(paste0(path, "D1112.XPT")) %>% mutate(cycle = "2011-2012"),
  by = c("SEQN", "cycle")
)

data_1314 <- inner_join(
  read_xpt(paste0(path, "1314.XPT")) %>% mutate(cycle = "2013-2014"),
  read_xpt(paste0(path, "D1314.XPT")) %>% mutate(cycle = "2013-2014"),
  by = c("SEQN", "cycle")
)

combined_1114 <- bind_rows(data_1112, data_1314)

filtered_1114 <- combined_1114 %>%
  filter(RIAGENDR == 2) %>%   
  select(SEQN, cycle, RIDAGEYR, RIAGENDR, IMQ040, WTINT2YR, SDMVPSU, SDMVSTRA)

filtered_1114 <- filtered_1114 %>%
  mutate(WTINT8YR = WTINT2YR / 2)

filtered_1114_valid <- filtered_1114 %>%
  filter(
    RIDAGEYR >= 9 & RIDAGEYR <= 59,
    IMQ040 %in% c(1, 2)  
  )

svy_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTINT8YR,
  nest = TRUE,
  data = filtered_1114_valid
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



# Title: ManyEcoEvo analysis

# Project: Out of sample analysis

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(quantreg)
library(here)

# read in the euc_specification_data_wide.csv
euc_os <- read_csv(here("data/euc_specification_data_wide.csv"))
glimpse(euc_os)

# calculate the three composite grass variables used:
# 1. exotic grass cover
# 2. native grass cover
# 3. total grass cover

# create an exotic grass cover variable
euc_os$exotic_grass <- 
  euc_os %>%
  select("ExoticAnnualGrass_cover", "ExoticPerennialGrass_cover") %>%
  rowSums()

# create a native grass variable
euc_os$native_perennial_grass <- 
  euc_os %>%
  select("NativePerennialGrass_cover", "NativePerennialGraminoid_cover") %>%
  rowSums()

# create a total grass variable
euc_os$total_grass <- 
  euc_os %>%
  select("NativePerennialGrass_cover", "NativePerennialGraminoid_cover",
         "ExoticAnnualGrass_cover", "ExoticPerennialGrass_cover"
  ) %>%
  rowSums()


# summarise these variables along with TAP at the property-scale
nrow(euc_os)
names(euc_os)

euc_os.m <- 
  euc_os %>%
  group_by(Property, Season) %>%
  summarise_at(vars(all_of( c("exotic_grass", "native_perennial_grass", "total_grass", "annual_precipitation") ) ),
               ~ mean(., na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Property) %>%
  summarise_at(vars(all_of( c("exotic_grass", "native_perennial_grass", "total_grass", "annual_precipitation") ) ), 
               list( ~ mean(., na.rm = TRUE),
                     ~ sd(., na.rm = TRUE) ) )


# summary of the models
summary(lm_out_glm[[3]])
lm_out_glm[[3]]$df.residual
summary(lm_out_glm[[4]])
lm_out_glm[[4]]$df.re
summary(lm_out_glm[[5]])
lm_out_glm[[5]]$df.re


# five models from Table 2
ci.out <- 
  lapply(lm_out_glm[2:5], function(model) {
  
  df.m <- predict(object = model, newdata = euc_os.m, se.fit = TRUE)
  
  data.frame(SurveyID = c("Q1", "Q2", "Q3"),
             fit = df.m$fit,
             se.fit = df.m$se.fit,
             ci.low = df.m$fit - (1.96*df.m$se.fit),
             ci.high = df.m$fit - (1.96*df.m$se.fit))
  
})

# write these to .csv files
for (i in 1:length(ci.out)) {
  
  write_csv(x = ci.out[[i]], 
            file = paste(here("figures_tables"), "/JamesHagan-predictions.", i, ".csv", sep = "") )
  
}

# get predictions from the two quantile regression models and output to .csv
predict(object = qm_1, newdata = euc_os.m, interval = "confidence", level = 0.95) %>%
  as.data.frame(.) %>%
  mutate(SurveyID = c("Q1", "Q2", "Q3")) %>%
  select(SurveyID, fit, ci.low = lower, ci.high = higher) %>%
  write_csv(x = ., file = here("figures_tables/quantile_regression_prediction_1.csv"))

# get t-value from coefficient that best answers the question
summary.rq(qm_1, se = "boot", bsmethod= "wild")
qm_1

predict(object = qm_2, newdata = euc_os.m, interval = "confidence", level = 0.95) %>%
  as.data.frame(.) %>%
  mutate(SurveyID = c("Q1", "Q2", "Q3")) %>%
  select(SurveyID, fit, ci.low = lower, ci.high = higher) %>%
  write_csv(x = ., file = here("figures_tables/quantile_regression_prediction_2.csv"))

# get t-value from coefficient that best answers the question
summary.rq(qm_2, se = "boot", bsmethod= "wild")




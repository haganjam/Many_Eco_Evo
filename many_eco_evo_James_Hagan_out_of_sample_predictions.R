
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


# five models from Table 2


y <- lm_out_glm[[2]]
y

df.m <- predict(object = lm_out_glm[[2]], newdata = euc_os.m, se.fit = TRUE)

data.frame(SurveyID = c("Q1", "Q2", "Q3"),
           fit = df.m$fit,
           se.fit = df.m$se.fit,
           ci.low = df.m$fit - (1.96*df.m$se.fit),
           ci.high = df.m$fit - (1.96*df.m$se.fit))






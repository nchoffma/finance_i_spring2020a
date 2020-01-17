
# Data and packages -------------------------------------------------------

library(tidyverse)
library(wrapr)
library(stringr)
library(lubridate)
inds <- read.csv("Industry49_data.csv", header = T)


# Question 5 --------------------------------------------------------------

# Clean the excess returns data
ers <- inds %>% 
  mutate(date = paste(year, month, day, sep = "-"),
         date = ymd(date)) %>% 
  select(date, starts_with("eR")) %>% 
  pivot_longer(-date, names_to = "industry", values_to = "eR") %>% 
  mutate(industry = str_remove_all(industry, "eR_"),
         industry = str_remove_all(industry, "I_"))

# a)
exp_ers <- ers %>% 
  group_by(industry) %>% 
  summarise(
    r_bar = mean(eR, na.rm = T),
    r_bar = r_bar * 12 + 0.01 # convert to ann, add RfR
  )

R_bar <- exp_ers[["r_bar"]] # get vector, which is what we're after
names(R_bar) <- exp_ers[["industry"]] # want to keep track

# b)
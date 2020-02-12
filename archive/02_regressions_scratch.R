
# Data and Packages -------------------------------------------------------

library(tidyverse)
library(wrapr)
library(stringr)
library(lubridate)
library(stargazer)
inds <- read.csv("Industry49_data.csv", header = T)


# Question 2 --------------------------------------------------------------

# Clean the data for regressions
ers <- inds %>% 
  mutate(date = paste(year, month, day, sep = "-"),
         date = ymd(date)) %>% 
  select(date, starts_with("eR"))

# Run multiple regression
y <- select(ers, -c(date, eR_Market)) %>% as.matrix()
rm <- ers["eR_Market"] %>% as.matrix()
models <- lm(formula = y ~ rm)

# Get coefficient values
coeffs <- data.frame(t(models$coefficients))
names(coeffs) <- qc(alpha, beta)
coeffs["industry"] <- rownames(coeffs)
rownames(coeffs) <- c()

# Get p-values
all_coefs <- coef(summary(models))
all_pvals <- data.frame(t(sapply(all_coefs, function(m) m[, 4])))
names(all_pvals) <- qc(p_alpha, p_beta)
all_pvals["industry"] <- str_remove_all(rownames(all_pvals), "Response ")
rownames(all_pvals) <- c()

# merge results
capm_results <- merge(coeffs, all_pvals, by = "industry")
capm_results <- mutate(capm_results, 
                       industry = str_remove_all(industry, "eR_I_"))

# export results
capm_results_out <- capm_results %>% 
  select(industry, beta, alpha, p_alpha) %>% 
  stargazer(., summary = F)

write(capm_results_out, file = "assignment_writeups/02_assign/table_1_capm.txt")

# One issue: do we have autocorrelation here? may need to try something 
# different (MLE, Sandwich, etc.) to get the right standard errors.
# The estimates of alpha should be consistent, though.

# B 
capm_results %>% 
  ggplot()

# Question 3 --------------------------------------------------------------

P <- c(1.575, 1.35, 3.425)
D <- matrix(c(2, 1, 4, 1, 3, 3), nrow = 2, byrow = T)


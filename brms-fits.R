library(cmdstanr)
library(brms)
library(tidybayes)
library(modeldb)

options(mc.cores = 4,
        brms.backend = "cmdstanr")

db <- d %>%
  add_dummy_variables(cohort_year, 
    values=c(-Inf,2014,2016), remove_original = F) %>%
  add_dummy_variables(year, 
    values=c(2011:2020), remove_original = F)

b1 <-
  brm(data = db, 
      family = gaussian(),
      y ~ 1 + (1 | state) + cohort_year_2014*year_2014 + 
               cohort_year_2014*year_2015 +
               cohort_year_2014*year_2016 +
               cohort_year_2014*year_2017 +
               cohort_year_2014*year_2018 +
               cohort_year_2014*year_2019 +
               cohort_year_2014*year_2020 +
               cohort_year_2016*year_2016 +
               cohort_year_2016*year_2017 +
               cohort_year_2016*year_2018 +
               cohort_year_2016*year_2019 +
               cohort_year_2016*year_2020 +
               cohort_year_2014 + cohort_year_2016 +
               year_2012 + year_2013 +
               year_2014 + year_2015 +
               year_2016 + year_2017 +
               year_2018 + year_2019 +
               year_2020,
      prior = c(prior(normal(0, 1.5), class = Intercept),  # bar alpha
                prior(normal(0, 0.5), class = b),          # betas
                prior(exponential(1), class = sd)),        # sigma
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes")

b1a <- update(b1,
         prior = c(prior(normal(0, 1.5), class = Intercept),
                  prior(normal(0, 0.5), class = b),                                         prior(exponential(1), class = sd)))


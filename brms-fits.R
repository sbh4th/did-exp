library(cmdstanr)
library(brms)
library(tidybayes)
library(modeldb)
library(ggplot2)

options(mc.cores = 4,
        brms.backend = "cmdstanr")

db <- d %>%
  add_dummy_variables(cohort_year, 
    values=c(-Inf,2014,2016), remove_original = F) %>%
  add_dummy_variables(year, 
    values=c(2011:2020), remove_original = F)

## Model without random effects at state level
b1 <-
  brm(data = db, 
      family = gaussian(),
      y ~ 1 + cohort_year_2014*year_2014 + 
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
      prior = c(prior(normal(0, 100), class = Intercept),
        prior(normal(0, 100), class = b)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes",
      file = "brms-fe")


## now model including random effects for state

b2 <-
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
      prior = c(prior(normal(0, 100), class = Intercept),
        prior(normal(0, 100), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes",
      file = "brms-re-flat")

b2_p <- prior_draws(b2) 

b2_p %>% ggplot(aes(x = Intercept)) + 
  stat_halfeye(color= "black", fill =  '#1b9e77') +
  labs(x = "Parameter value", 
    y = "Density") +
  theme_classic()

b2 %>% 
  mcmc_plot(variable = "^r_", regex = T) +
  theme(axis.text.y = element_text(hjust = 0))

b3 <-
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
      prior = c(prior(normal(0, 10), class = Intercept),
        prior(normal(0, 5), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes",
      file = "brms-re-inf")

# compare models using leave-one-out criterion
b1 <- add_criterion(b1, "loo")
b2 <- add_criterion(b2, "loo")
b3 <- add_criterion(b3, "loo")
b5 <- add_criterion(b5, "loo")

l <- loo_compare(b1, b2, b3, b5, criterion = "loo")
print(l, simplify = F)

## all treated (treated = 1) and all untreated (treated = 0)
db_0 <- db %>% 
  mutate(treat = 0)
db_1 <- db %>% 
  mutate(treat = 1)
db_me <- bind_rows(db_0, db_1)

# marginal predicted risks
mpy <- add_predicted_draws(b3, newdata=db_me,
    allow_new_levels = TRUE, ndraws=10) %>%
  mutate(Tx = recode_factor(treat, 
    `0` = "Control", `1` = "Treated")) %>%
  group_by(Tx, .draw) %>%
  summarise(`Pr(y)` = mean(`.prediction`))




modelsummary(list("Fixed only" = b1, 
 "State RE (non-informative)" = b2, 
 "State RE (informative)" = b3), 
  coef_map = cm, statistic="mad",
  shape = term ~ model,
  gof_omit ='._*')




b5 <-
  brm(data = db, 
      family = gaussian(),
      y ~ 1 + (1 | state) + treat*cohort_year_2014*year_2014 + 
               treat*cohort_year_2014*year_2015 +
               treat*cohort_year_2014*year_2016 +
               treat*cohort_year_2014*year_2017 +
               treat*cohort_year_2014*year_2018 +
               treat*cohort_year_2014*year_2019 +
               treat*cohort_year_2014*year_2020 +
               treat*cohort_year_2016*year_2016 +
               treat*cohort_year_2016*year_2017 +
               treat*cohort_year_2016*year_2018 +
               treat*cohort_year_2016*year_2019 +
               treat*cohort_year_2016*year_2020 +
               cohort_year_2014 + cohort_year_2016,
      prior = c(prior(normal(0, 10), class = Intercept),
        prior(normal(0, 5), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4)

b6 <-
  brm(data = db, 
      family = gaussian(),
      y ~ 1 + (1 | state) + treat:cohort_year_2014:year_2014 +
        cohort_year_2014 + year_2014,
      iter = 2000, warmup = 1000, chains = 4, cores = 4)

b7 <-
  brm(data = db, 
      family = gaussian(),
      y ~ 1 + (1 | state) + treat:cohort_year_2014:year_2014 + 
               treat:cohort_year_2014:year_2015 +
               treat:cohort_year_2014:year_2016 +
               treat:cohort_year_2014:year_2017 +
               treat:cohort_year_2014:year_2018 +
               treat:cohort_year_2014:year_2019 +
               treat:cohort_year_2014:year_2020 +
               treat:cohort_year_2016:year_2016 +
               treat:cohort_year_2016:year_2017 +
               treat:cohort_year_2016:year_2018 +
               treat:cohort_year_2016:year_2019 +
               treat:cohort_year_2016:year_2020 +
               cohort_year_2014 + cohort_year_2016 +
               year_2012 + year_2013 +
               year_2014 + year_2015 +
               year_2016 + year_2017 +
               year_2018 + year_2019 +
               year_2020,
      prior = c(prior(normal(0, 10), class = Intercept),
        prior(normal(0, 5), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4)

bme_c <- slopes(
  b5, 
  newdata   = subset(db, treat & cohort_year),
  variables = "treat", 
  by        = "cohort_year"
  )

bme_t <- slopes(
  b5, 
  newdata   = subset(db, treat & year>=2014),
  variables = "treat", 
  by        = "year"
  )

d <- d %>%
  mutate(event = year - cohort_year)

# Group effects for ETWFE (by time since treatment)
me_e <- slopes(
  etwfe, 
  newdata   = subset(d, treat & event>=0),
  variables = "treat", 
  by        = "event"
  )

# marginal effects 

bme_avg <- slopes(
  b5, 
  newdata   = subset(db, treat==1),
  variables = "treat", 
  by        = "treat"
  )

# wrangle aggregate ATTs for model summary table
bti <- data.frame(
  term = paste("ATT(", bme_avg$term, ")", sep = ""),
  estimate = bme_avg$estimate,
  conf.low = bme_avg$conf.low,
  conf.high = bme_avg$conf.high,
  std.error = abs(bme_avg$conf.high - 
                    bme_avg$conf.low) / (2 * 1.96)
)

gl <- data.frame()

etwfe_me_avg <- list(tidy = ti, glance = gl)
class(etwfe_me_avg) <- "modelsummary_list"

modelsummary(list("Simple Average" = etwfe_me_avg),
  shape = term ~ model + statistic, 
  statistic = c("({std.error})", "conf.int"),
  gof_omit ='._*')


# rename coefficients from the models for table
cm5 <- c("b_cohort_year_2014:year_2014" = "ATT(2014,2014)",
      "b_treat:cohort_year_2014:year_2015" = "ATT(2014,2015)",
      "b_treat:cohort_year_2014:year_2016" = "ATT(2014,2016)",
      "b_treat:cohort_year_2014:year_2017" = "ATT(2014,2017)",
      "b_treat:cohort_year_2014:year_2018" = "ATT(2014,2018)",
      "b_treat:cohort_year_2014:year_2019" = "ATT(2014,2019)",
      "b_treat:cohort_year_2014:year_2020" = "ATT(2014,2020)",
      "b_treat:year_2016:cohort_year_2016" = "ATT(2016,2016)",
      "b_treat:year_2017:cohort_year_2016" = "ATT(2016,2017)",
      "b_treat:year_2018:cohort_year_2016" = "ATT(2016,2018)",
      "b_treat:year_2019:cohort_year_2016" = "ATT(2016,2019)",
      "b_treat:year_2020:cohort_year_2016" = "ATT(2016,2020)")

# build the table
modelsummary(b5, 
  coef_map = cm5, statistic="mad",
  gof_omit ='._*')

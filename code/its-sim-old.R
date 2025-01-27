# code for experimenting with ITS simulation. 
library(here)
library(tidyverse)
library(estimatr)
library(modelsummary)

set.seed(48620)

# unit fixed effects (unobserved heterogeneity)
unit <- tibble(
  unit = 1:1000,
  # generate clusters
  state = sample(1:30, 1000, replace = TRUE),
  unit_fe = rnorm(1000, state/10, 1),
  # generate instantaneous treatment effect
  #mu = rnorm(nobs, true_mu, 0.2)
  mu = 2
)

# year fixed effects (first part)
year <- tibble(
  year = 2001:2020,
  year_fe = rnorm(length(year), 0, 0.5)
)

# Put the clusters into treatment groups
treat_taus <- tibble(
  # sample the clusters randomly
  state = sample(1:30, 30, replace = FALSE),
  # place the randomly sampled states into 1\{t \ge g \}G_g
  cohort_year = sort(rep(c(2010, 2015), times=c(10,20)))
)

# make main dataset
# full interaction of unit X year 
dtest <- expand_grid(unit = 1:1000, year = 2001:2020) %>% 
  left_join(., unit) %>% 
  left_join(., year) %>% 
  left_join(., treat_taus) %>% 
  # make error term and get treatment indicators and treatment effects
  # Also get cohort specific trends (modify time FE)
  mutate(error = rnorm(1000*20, 0, 2),
         treat = ifelse((year >= 
           cohort_year), 1, 0),
         # treatment effect = 3 if 2015, 6 if 2010, annually
         mu = ifelse(cohort_year==2015, -2, 2),
         tau = ifelse(treat == 1, mu, 0),
         # year trends differ by cohort
         year_fe = ifelse(cohort_year==2015 & year<2005,
          year_fe + (2015 - cohort_year) * 
           (year - cohort_year) / 5 + 12,
          year_fe + (2020 - cohort_year) *
           (year - cohort_year) / 5 + 20)
  ) %>% 
  # calculate cumulative treatment effects
  group_by(unit) %>% 
  mutate(tau_cum = cumsum(tau)) %>% 
  ungroup() %>% 
  # calculate the dependent variable
  mutate(y = (2020 - cohort_year) + 
           unit_fe + year_fe + tau_cum + error) 
  # Relabel 2018 cohort as never-treated
  # mutate(cohort_year = ifelse(cohort_year == 2018, Inf, cohort_year))

# bring in simulated DiD dataset
dss <- dtest %>%
  # read_rds(here("data", 
 #"did-sim-data.rds")) %>%
  mutate(year0 = year - 2001,
    yearc = year - cohort_year,
    cohort_2015 = ifelse(cohort_year == 2015, 1, 0))

# plots
dss %>% 
  ggplot(aes(x = year0, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = dss %>% 
    group_by(cohort_year, year0) %>% 
    summarize(y = mean(y)),
   aes(x = year0, y = y, group = factor(cohort_year),
    color = factor(cohort_year)), linewidth = 2) + 
  geom_line(data = dss %>%
    group_by(year0) %>%
      summarize(y = mean(y)),
    aes(x = year0, y = y, color = "#4daf4a", group="All"), 
    linewidth = 1, linetype = 2) +
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(5, 10, 15, 20)) +
  geom_vline(xintercept = 10, color = '#E41A1C', linewidth = 2) + 
  geom_vline(xintercept = 15, color = '#377EB8', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1', 
    labels = c("2010", "2015", "Average")) + 
  theme_classic()

# plot of average (ignoring cohorts)
dss %>% 
  ggplot(aes(x = year0, y = y)) + 
  geom_line(aes(group=unit), alpha = 1/8, color = "grey") +
  geom_line(data = ds %>% 
    group_by(year0) %>% 
    summarize(y = mean(y)),
   aes(x = year0, y = y), linewidth=2) + 
  labs(x = "", y = "Y") +
  scale_x_continuous(breaks = c(5, 10, 15, 20)) 

# relative time
dss %>% 
  ggplot(aes(x = yearc, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = dss %>% 
    group_by(cohort_year, yearc) %>% 
    summarize(y = mean(y)),
   aes(x = yearc, y = y, group = factor(cohort_year),
    color = factor(cohort_year)), linewidth = 2) + 
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(-5, 0, 5, 10)) +
  geom_vline(xintercept = 0, color = '#E41A1C', linewidth = 2) + 
  geom_vline(xintercept = 0, color = '#377EB8', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1') 

# relative time
dss %>% 
  ggplot(aes(x = yearc, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = dss %>% 
    group_by(cohort_year, yearc) %>% 
    summarize(y = mean(y)),
   aes(x = yearc, y = y, group = factor(cohort_year),
    color = factor(cohort_year)), linewidth = 2) +
  geom_line(data = dss %>%
    group_by(yearc) %>%
      summarize(y = mean(y)),
    aes(x = yearc, y = y, color = "#4daf4a", group="All"), 
    linewidth = 2) +
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(-5, 0, 5, 10)) +
  geom_vline(xintercept = 0, color = '#984ea3', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1', labels = c("2010", "2015", "Average")) + theme_classic()

# relative time
dss %>%
  filter(yearc >= -9 & yearc <= 5) %>%
  ggplot(aes(x = yearc, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = dss %>% 
    filter(yearc >= -9 & yearc <= 5) %>%
    group_by(cohort_year, yearc) %>% 
    summarize(y = mean(y)),
   aes(x = yearc, y = y, group = factor(cohort_year),
    color = factor(cohort_year)), linewidth = 2) +
  geom_line(data = dss %>%
    filter(yearc >= -9 & yearc <= 5) %>%
    group_by(yearc) %>%
      summarize(y = mean(y)),
    aes(x = yearc, y = y, color = "#4daf4a", group="All"), 
    linewidth = 2) +
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(-5, 0, 5, 10)) +
  geom_vline(xintercept = 0, color = '#984ea3', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1', labels = c("2010", "2015", "Average")) + theme_classic()

m1 <- lm_robust(y ~ year, data=dtest, clusters = state)
tidy(m1)

m2 <- lm_robust(y ~ year + cohort_year, 
  data=dtest, clusters = state)
tidy(m2)

m3 <- lm_robust(y ~ year + cohort_year + 
  year:cohort_year, data=dtest, clusters = state)
modelsummary(list(m1,m2,m3))

mb <- lm_robust(y ~ year0 + treat,
  data = ds, clusters = state)
mbs <- lm_robust(y ~ year0 + treat + factor(cohort_year),
  data = ds, clusters = state)
mbc <- lm_robust(y ~ yearc + treat,
  data = ds, clusters = state)
mbcs <- lm_robust(y ~ yearc + treat + factor(cohort_year),
  data = ds, clusters = state)
mbcsr <- lm_robust(y ~ yearc + treat + factor(cohort_year),
  data = subset(ds, yearc>=-9 & yearc<=5), 
  clusters = state)

# ITS model with slope and level change
mits <- lm_robust(y ~ year0 + treat + 
  yearct, data = ds, clusters = state)
mitsc <- lm_robust(y ~ year0 + treat + 
  yearct + factor(cohort_year), 
  data = ds, clusters = state)

# ITS model with centered year
mits_yc <- lm_robust(y ~ yearc + treat + 
  yearct, data = ds, clusters = state)
mitsc_yc <- lm_robust(y ~ year0 + treat + 
  yearct + factor(cohort_year), 
  data = ds, clusters = state)

# ITS model with centered year, restricted
mits_yc_r <- lm_robust(y ~ yearc + treat + 
  yearct, data = subset(ds, yearc>=-9 & yearc<=5), 
  clusters = state)
mitsc_yc_r <- lm_robust(y ~ yearc + treat + 
  yearct + factor(cohort_year), 
  data = subset(ds, yearc>=-9 & yearc<=5),
  clusters = state)

modelsummary(list(mits,mitsc))


ds %>% 
  filter(cohort_2015==1) %>%
  ggplot(aes(x = year0, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_smooth(data = ds %>%
    filter(cohort_2015==1) %>%
    group_by(year0, treat) %>%
      summarize(y = mean(y)),
    aes(x = year0, y = y, color = factor(treat), 
        group = factor(treat)), 
    linewidth = 1, method='lm') +
  geom_smooth(data = ds %>%
    filter(cohort_2015==1) %>%
    group_by(year0, treat) %>%
      summarize(y = mean(y)),
    aes(x = year0, y = y, color = factor(treat), 
        group = factor(treat)), 
    linewidth = 1, method='lm') +
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(5, 10, 15, 20)) +
  geom_vline(xintercept = 15, color = '#377EB8', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1', 
    labels = c("Pre", "Post")) + 
  theme_classic()

dss %>% 
  filter(cohort_2015==1) %>%
  ggplot(aes(x = year0, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_smooth(data = dss %>%
    filter(cohort_2015==1) %>%
    group_by(year0, treat) %>%
      summarize(y = mean(y)),
    aes(x = year0, y = y, color = factor(treat), 
        group = factor(treat)), 
    linewidth = 1, method='lm') +
  geom_smooth(data = dss %>%
    filter(cohort_2015==1) %>%
    group_by(year0, treat) %>%
      summarize(y = mean(y)),
    aes(x = year0, y = y, color = factor(treat), 
        group = factor(treat)), 
    linewidth = 1, method='lm') +
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(5, 10, 15, 20)) +
  geom_vline(xintercept = 15, color = '#377EB8', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1', 
    labels = c("Pre", "Post")) + 
  theme_classic()

# experimental code
eof_predictions <- function(model, newdata) {
  
  # Generate predictions
  pred <- predictions(model, newdata = newdata)
  
  # Compute the hypothesis test for the difference
  hypothesis_test <- predictions(model, newdata = newdata, hypothesis = "b2 - b1 = 0")
  
  # Return both results as a list
  list(predictions = pred, hypothesis_test = hypothesis_test)
}


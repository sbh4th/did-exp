# load packages
library(here)
library(tidyverse)
library(estimatr)
library(modelsummary)
library(marginaleffects)
library(tinytable)

# set seed for reproducibility
set.seed(95478)

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
         # treatment effect = 4 if 2015, 1 if 2010, annually
         mu = ifelse(cohort_year==2015, 4, 1),
         tau = ifelse(treat == 1, mu, 0),
         # year trends differ by cohort
         year_fe = year_fe + (2020 - cohort_year) * 
           (year - cohort_year) / 5 + 20
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
ds <- dtest %>%
  # read_rds(here("data", 
 #"did-sim-data.rds")) %>%
  mutate(year0 = year - 2001,
    yearc = year - cohort_year,
    since = treat * yearc,
    cohort_2015 = ifelse(cohort_year == 2015, 1, 0))

ds %>% 
  ggplot(aes(x = year0, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = ds %>% 
    group_by(cohort_year, year0) %>% 
    summarize(y = mean(y)),
   aes(x = year0, y = y, group = factor(cohort_year),
    color = factor(cohort_year)), linewidth = 2) + 
  geom_line(data = ds %>%
    group_by(year0) %>%
      summarize(y = mean(y)),
    aes(x = year0, y = y, color = "#4daf4a", group="All"), 
    linewidth = 1, linetype = 2) +
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(5, 10, 15, 20)) +
  geom_vline(xintercept = 9, color = '#E41A1C', linewidth = 2) + 
  geom_vline(xintercept = 14, color = '#377EB8', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1', 
    labels = c("2010", "2015", "Average")) + 
  theme_classic()

# models
m_c2010 <- lm_robust(y ~ year0 + treat + since, 
  data = subset(ds, cohort_year == 2010), 
  clusters = state)

m_c2015 <- lm_robust(y ~ year0 + treat + since, 
  data = subset(ds, cohort_year == 2015), 
  clusters = state)

modelsummary(list("2010 cohort" = m_c2010,
                  "2015 cohort" = m_c2015),
  fmt = 2, statistic = "conf.int", 
  shape = term ~ model + statistic,
  gof_omit = 'DF|Deviance|R2|AIC|BIC|RMSE')

# predictions for 2010 cohort
p_2010 <- predictions(m_c2010, 
  newdata = datagrid(year0 = c(19,19),
    treat = c(0,1), since = c(0,10))) 

d_2010 <- hypotheses(p_2010, "b4 - b1 = 0")

p_2010 <- p_2010 %>%
  filter(rowid==c(1,4)) %>%
  bind_rows(d_2010) %>%
  mutate(term = c("Untreated", "Treated",
    "Difference"))

# predictions for 2015 cohort
p_2015 <- predictions(m_c2015, 
  newdata = datagrid(year0 = c(19,19),
    treat = c(0,1), since = c(0,5)))
d_2015 <- hypotheses(p_2015, "b4 - b1 = 0")

p_2015 <- p_2015 %>% 
    filter(rowid==c(1,4)) %>%
  bind_rows(d_2015) %>%
  mutate(term = c("Untreated", "Treated",
    "Difference"))

modelsummary(list("2010 cohort" = p_2010,
  "2015 cohort" = p_2015), fmt=2,
  statistic = "conf.int", shape = term ~ statistic)

te <- p_2010$estimate[3]*(1/3) + p_2015$estimate[3]*(2/3)


m_all <- lm_robust(formula = 
  y ~ year0 + treat + since, 
  data = ds, clusters = state)

m_all_s <- lm_robust(formula = 
  y ~ cohort_2015 + year0 + treat + since, 
  data = ds, clusters = state)

m_all_int <- lm_robust(formula = 
  y ~ cohort_2015 * (year0 + treat + since), 
  data = ds, clusters = state)

modelsummary(list("2010 cohort" = m_c2010,
  "2015 cohort" = m_c2015, "No cohort FE" = m_all,
  "Cohort FE" = m_all_s, "Interactive" = m_all_int),
  fmt = 2, gof_omit = 'DF|Deviance|R2|AIC|BIC|RMSE')

# predictions at end of follow-up by treatment
# set up dataframe for predictions
eof <- tibble(
  cohort_2015 = c(0,0,1,1),
  year0 = c(rep(19,4)),
  treat = c(0,1,0,1),
  since = c(0,10,0,5)
)

m_all_int_ap <- predictions(m_all_int, 
  newdata = eof)

# treatment effects
m_all_int_te <- hypotheses(m_all_int_ap, 
  c("b2 - b1 = 0", "b4 - b3 = 0")) %>%
  mutate(term = c("2010 cohort", "2015 cohort"))

modelsummary(list("Marginal effects from interaction model" 
  = m_all_int_te), fmt=2,
  statistic = "conf.int", shape = term ~ statistic,
  notes = "Note: calculated at end of follow-up.")


m_all_int_slopes <- avg_slopes(m_all_int, 
  variable=c("treat","since"), 
  by="cohort_2015") %>%
  mutate(
    term = c("Post-slope (since)","Post-slope (since)",
             "Level shift (treat)","Level shift (treat)"),
    cohort_2015 = factor(cohort_2015,
    labels = c("2010 cohort","2015 cohort")))

modelsummary(list("Average slopes: Interaction model" 
                  = m_all_int_slopes),
  fmt = 2, statistic = "conf.int", 
  shape = term ~ cohort_2015 + statistic,
  gof_omit = 'DF|Deviance|R2|AIC|BIC|RMSE')

# predictions at end of follow-up
# averaging over the entire sample
m_all_int_aap <- hypotheses(m_all_int_ap, 
  c("(1/3)*b1 + (2/3)*b3 = 0", 
  "(1/3)*b2 + (2/3)*b4 = 0")) %>%
  mutate(term = c("Untreated", "Treated"))

m_all_int_ate <- hypotheses(m_all_int_ap, 
  "(1/3)*(b2 - b1) + (2/3)*(b4 - b3) = 0") %>%
  mutate(term = "Difference")

m_all_int_ates <- m_all_int_aap %>% 
  bind_rows(m_all_int_ate)

modelsummary(list("Population average effect: interaction model" 
    = m_all_int_ates),
  fmt = 2, statistic = "conf.int", 
  shape = term ~ model + statistic,
  gof_omit = 'DF|Deviance|R2|AIC|BIC|RMSE',
  notes = "Note: cohort fixed effect set to 2/3.")


ds %>%
  ggplot(aes(x = yearc, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = ds %>% 
    group_by(cohort_year, yearc) %>% 
    summarize(y = mean(y)),
   aes(x = yearc, y = y, group = factor(cohort_year),
    color = factor(cohort_year)), linewidth = 2) +
  geom_line(data = ds %>%
    group_by(yearc) %>%
      summarize(y = mean(y)),
    aes(x = yearc, y = y, color = "#4daf4a", group="All"), 
    linewidth = 1, linetype = 2) +
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(-5, 0, 5, 10)) +
  geom_vline(xintercept = 0, color = '#984ea3', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1', labels = c("2010", "2015", "Average")) + theme_classic()

m_c2010_c <- lm_robust(
  y ~ yearc + treat + since, 
  data = subset(ds, cohort_2015==0), 
  clusters = state)

m_c2015_c <- lm_robust(
  y ~ yearc + treat + since, 
  data = subset(ds, cohort_2015==1), 
  clusters = state)

m_all_c <- lm_robust(
  y ~ yearc + treat + since, 
  data = ds, clusters = state)

m_all_s_c <- lm_robust(
  y ~ cohort_2015 + yearc + treat + since, 
  data = ds, clusters = state)

m_all_int_c <- lm_robust(
  y ~ cohort_2015 * (yearc + treat + since), 
  data = ds, clusters = state)

modelsummary(list("2010 cohort" = m_c2010_c,
  "2015 cohort" = m_c2015_c, "No cohort FE" = m_all_c,
  "Cohort FE" = m_all_s_c, "Interactive" = m_all_int_c),
  fmt = 2, gof_omit = 'DF|Deviance|R2|AIC|BIC|RMSE')

# predictions for 2010 cohort
p_2010_c <- predictions(m_c2010_c, 
  newdata = datagrid(yearc = c(10,10),
    treat = c(0,1), since = c(0,10))) 

d_2010_c <- hypotheses(p_2010_c, "b4 - b1 = 0")

p_2010_c <- p_2010_c %>%
  filter(rowid==c(1,4)) %>%
  bind_rows(d_2010_c) %>%
  mutate(term = c("Untreated", "Treated",
    "Difference"))

# predictions for 2015 cohort
p_2015_c <- predictions(m_c2015_c, 
  newdata = datagrid(yearc = c(5,5),
    treat = c(0,1), since = c(0,5)))
d_2015_c <- hypotheses(p_2015_c, "b4 - b1 = 0")

p_2015_c <- p_2015_c %>% 
    filter(rowid==c(1,4)) %>%
  bind_rows(d_2015_c) %>%
  mutate(term = c("Untreated", "Treated",
    "Difference"))

modelsummary(list("2010 cohort (at 10yrs)" = p_2010_c,
  "2015 cohort (at 5 yrs)" = p_2015_c), fmt=2,
  statistic = "conf.int", shape = term ~ statistic)

# relative time with restriction
ds %>%
  filter(yearc >= -9 & yearc <= 5) %>%
  ggplot(aes(x = yearc, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = ds %>% 
    filter(yearc >= -9 & yearc <= 5) %>%
    group_by(cohort_year, yearc) %>% 
    summarize(y = mean(y)),
   aes(x = yearc, y = y, group = factor(cohort_year),
    color = factor(cohort_year)), linewidth = 2) +
  geom_line(data = ds %>%
    filter(yearc >= -9 & yearc <= 5) %>%
    group_by(yearc) %>%
      summarize(y = mean(y)),
    aes(x = yearc, y = y, color = "#4daf4a", group="All"), 
    linewidth = 1, linetype = 2) +
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(-5, 0, 5, 10)) +
  geom_vline(xintercept = 0, color = '#984ea3', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1', labels = c("2010", "2015", "Average")) + theme_classic()

m_c2010_cr <- lm_robust(
  y ~ yearc + treat + since, 
  data = subset(ds, cohort_2015==0 &
    yearc >= -9 & yearc <= 5), 
  clusters = state)

m_c2015_cr <- lm_robust(
  y ~ yearc + treat + since, 
  data = subset(ds, cohort_2015==1 &
    yearc >= -9 & yearc <= 5), 
  clusters = state)

m_all_cr <- lm_robust(
  y ~ yearc + treat + since, 
  data = subset(ds, yearc >= -9 & yearc <= 5), 
  clusters = state)

m_all_s_cr <- lm_robust(
  y ~ cohort_2015 + yearc + treat + since, 
  data = subset(ds, yearc >= -9 & yearc <= 5), 
  clusters = state)

m_all_int_cr <- lm_robust(
  y ~ cohort_2015 * (yearc + treat + since), 
  data = subset(ds, yearc >= -9 & yearc <= 5), 
  clusters = state)

modelsummary(list("2010 cohort" = m_c2010_cr,
  "2015 cohort" = m_c2015_cr, "No cohort FE" = m_all_cr,
  "Cohort FE" = m_all_s_cr, "Interactive" = m_all_int_cr),
  fmt = 2, gof_omit = 'DF|Deviance|R2|AIC|BIC|RMSE')

ds2 <- ds %>% 
  # add non-linearity to early part of 2015 cohort
  mutate(
    year_fe = ifelse(cohort_year==2015 & year<2006,
          year_fe + (2015 - cohort_year) * 
           (year - cohort_year) / 5 + 12,
          year_fe + (2020 - cohort_year) *
           (year - cohort_year) / 5 + 20), 
    # re-generate outcome
    y = (2020 - cohort_year) + 
           unit_fe + year_fe + tau_cum + error) 

ds2 %>% 
  ggplot(aes(x = year0, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = ds2 %>% 
    group_by(cohort_year, year0) %>% 
    summarize(y = mean(y)),
   aes(x = year0, y = y, group = factor(cohort_year),
    color = factor(cohort_year)), linewidth = 2) + 
  geom_line(data = ds2 %>%
    group_by(year0) %>%
      summarize(y = mean(y)),
    aes(x = year0, y = y, color = "#4daf4a", group="All"), 
    linewidth = 1, linetype = 2) +
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(5, 10, 15, 20)) +
  geom_vline(xintercept = 9, color = '#E41A1C', linewidth = 2) + 
  geom_vline(xintercept = 14, color = '#377EB8', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1', 
    labels = c("2010", "2015", "Average")) + 
  theme_classic()

ds2 %>%
  filter(yearc >= -9 & yearc <= 5) %>%
  ggplot(aes(x = yearc, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = ds2 %>% 
    filter(yearc >= -9 & yearc <= 5) %>%
    group_by(cohort_year, yearc) %>% 
    summarize(y = mean(y)),
   aes(x = yearc, y = y, group = factor(cohort_year),
    color = factor(cohort_year)), linewidth = 2) +
  geom_line(data = ds2 %>%
    filter(yearc >= -9 & yearc <= 5) %>%
    group_by(yearc) %>%
      summarize(y = mean(y)),
    aes(x = yearc, y = y, color = "#4daf4a", group="All"), 
    linewidth = 1, linetype = 2) +
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(-5, 0, 5, 10)) +
  geom_vline(xintercept = 0, color = '#984ea3', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1', labels = c("2010", "2015", "Average")) + theme_classic()

# calendar time full model
m_all_int2 <- lm_robust(
  y ~ cohort_2015 * (year0 + treat + since),
  data = ds2, clusters = state)

# relative time full model
m_all_int_cr2 <- lm_robust(
  y ~ cohort_2015 * (yearc + treat + since), 
  data = subset(ds2, yearc >= -9 & yearc <= 5), 
  clusters = state)

ft <- modelsummary(list("Sim1: Interaction" = m_all_int, 
  "Sim2: Interaction" = m_all_int2,
  "Sim1: Interaction" = m_all_int_cr, 
  "Sim2: Interaction" = m_all_int_cr2),
  fmt = 2, coef_rename = c("year0" = "year", 
                  "yearc" = "year"),
  gof_omit = 'DF|Deviance|R2|AIC|BIC|RMSE',
  output = 'tinytable')

tt(ft@table_dataframe) %>%
  group_tt(j=list("Calendar time" = 2:3, 
    "Relative time (restricted)" = 4:5))

# relative time basic model
m_all_cr2 <- lm_robust(
  y ~ yearc + treat + since,
  data=subset(ds2, yearc >= -9 & yearc <= 5), 
  clusters = state)

# relative time cohort FE model
m_all_s_cr2 <- lm_robust(
  y ~ cohort_2015 + yearc + treat + since,
  data=subset(ds2, yearc >= -9 & yearc <= 5), 
  clusters = state)

modelsummary(list("No cohort FE" = m_all_cr2,
  "Cohort FE" = m_all_s_cr2),
  fmt = 2, statistic = "conf.int", 
  shape = term ~ model + statistic,
  gof_omit = 'DF|Deviance|R2|AIC|BIC|RMSE')


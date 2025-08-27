library(here)
library(tidyverse)
library(truncnorm)
library(etwfe)

# unit fixed effects (unobserved heterogeneity)
unit <- tibble(
  unit = 1:57000,
  # generate clusters
  country = sample(1:38, 57000, replace = TRUE),
  unit_fe = rnorm(57000, country/10, 1),
  # generate instantaneous treatment effect
  #mu = rnorm(nobs, true_mu, 0.2)
  # mu = 2
)

# wave fixed effects (first part)
wave <- tibble(
  wave = 1:4,
  wave_fe = rnorm(length(wave), 0, 0.5)
)

# Put the clusters into treatment groups
treat_taus <- tibble(
  # sample the clusters randomly
  country = sample(1:38, 38, replace = FALSE),
  # place the randomly sampled countries into 1\{t \ge g \}G_g
  cohort_wave = sort(c(rep(1, 24), rep(2, 2), 
    rep(3, 4), rep(4, 8)))
)

# make main dataset
# full interaction of unit X year 
d <- expand_grid(unit = 1:57000, wave = 1:4) %>% 
  left_join(., unit) %>% 
  left_join(., wave) %>% 
  left_join(., treat_taus) %>% 
  # make error term and get treatment indicators and treatment effects
  # Also get cohort specific trends (modify time FE)
  mutate(error = rnorm(57000*4, 0, 1),
    treat = ifelse((wave >= 
      cohort_wave)*(cohort_wave != 1), 1, 0),
      # treatment effect = 1 if 2016, 2 if 2014, annually
      mu = ifelse(cohort_wave==2, 0, 
        ifelse(cohort_wave==3, 0, 0)),
      tau = ifelse(treat == 1, mu, 0),
      wave_fe = wave_fe - 1.5*(wave - cohort_wave)
  ) %>% 
  # calculate cumulative treatment effects
  group_by(unit) %>% 
  mutate(tau_cum = cumsum(tau)) %>% 
  ungroup() %>% 
  # calculate the dependent variable
  mutate(y = rtruncnorm(length(unit),
    a = 0, b = 20,
    mean = 15 - cohort_wave + 
     unit_fe + wave_fe + tau_cum + error,
    sd = 4),
    y = round(y)) %>%
  # Relabel 2018 cohort as never-treated
  mutate(cohort_wave = ifelse(
    cohort_wave == 1, Inf, cohort_wave))
  
theme_set(theme_classic() + 
            theme(plot.background = element_blank()))

plot <- d %>% 
  ggplot(aes(x = wave, y = y, group = unit)) + 
  geom_line(alpha = 1/8, color = "grey") + 
  geom_line(data = d %>% 
    group_by(cohort_wave, wave) %>% 
    summarize(y = mean(y)),
   aes(x = wave, y = y, group = factor(cohort_wave),
    color = factor(cohort_wave)), linewidth = 2) + 
  labs(x = "", y = "Y",  color = "Treatment group   ") +
  scale_x_continuous(breaks = c(1, 2, 3, 4)) +
  geom_vline(xintercept = 2, 
    color = '#E41A1C', linewidth = 2) + 
  geom_vline(xintercept = 3, 
    color = '#377EB8', linewidth = 2) + 
  scale_color_brewer(palette = 'Set1') + 
  theme(legend.position = 'bottom',
        #legend.title = element_blank(), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_color_manual(labels = 
    c("2", "3", "4", "Never-treated"),
  values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984ea3")) +
  ggtitle("Simulated data with heterogeneous treatment effect dynamics across cohorts \n and with a never-treated group")+
  theme(plot.title = element_text(hjust = 0.5, size=12))

plot 

# fit the model
etwfe_model = fixest::feols(
  y ~ treat:i(cohort_wave, i.wave, ref=Inf, ref2 = 1) | 
    cohort_wave + wave,
  data = d,
  vcov = ~country
)

# averge marginal effect
etwfe_me <- slopes(etwfe_c, newdata = subset(d, treat==1), var = "treat", by = "treat")
  
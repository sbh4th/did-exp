library(tidyverse)
library(fixest)
library(patchwork)

set.seed(42)
N <- 2000

# ============================================================
# Setup / helpers
#
# Two-period differential exposure design:
#   Y_it = alpha_i + lambda_t + ATT(D_i) * Post_t + eps_it
#
# Parallel trends holds *exactly* by construction (lambda_t
# is common to all units; unit FEs are random noise).
#
# TWFE regression: feols(y ~ I(D * period) | unit + period)
#
# In the 2-period case, the TWFE coefficient reduces to:
#   beta_twfe = Cov(D_i, Delta_Y_i) / Var(D_i)
#             = slope of OLS regression of DeltaY on D
#
# This is NOT the average ATT = E[ATT(D_i)].
# Weights on unit i are proportional to (D_i - D_bar),
# so below-average-dose units get NEGATIVE weight.
# ============================================================

make_panel <- function(D, att_fn, time_fe = 1, sigma = 2) {
  n     <- length(D)
  tau   <- att_fn(D)
  alpha <- rnorm(n, mean = 50, sd = 10)
  # compute y outside tibble() to avoid tau column shadowing the local variable
  y_pre  <- alpha                      + rnorm(n, 0, sigma)
  y_post <- alpha + time_fe + tau      + rnorm(n, 0, sigma)
  tibble(
    unit   = rep(seq_len(n), 2),
    period = rep(0:1, each = n),
    D      = rep(D, 2),
    tau    = rep(tau, 2),
    y      = c(y_pre, y_post)
  )
}

twfe_coef <- function(panel) {
  fit <- feols(y ~ I(D * period) | unit + period,
               data = panel, warn = FALSE)
  unname(coef(fit)[1])
}

# ============================================================
# Scenario 1: Linear ATT
#   D ~ Uniform(5, 20)   [realistic ambient PM2.5 range]
#   ATT(d) = 0.5 * d     [each ug/m3 increases outcome by 0.5]
#
#   True avg ATT = 0.5 * E[D] = 0.5 * 12.5 = 6.25
#   TWFE estimate:  Cov(D, 0.5*D) / Var(D) = 0.5
#
#   TWFE gives the slope (a "causal response"), not the
#   average level effect. Different estimands — TWFE is
#   off by a factor of E[D] = 12.5.
# ============================================================
D1   <- runif(N, 5, 20)
att1 <- function(d) 0.5 * d

p1       <- make_panel(D1, att1)
beta1    <- twfe_coef(p1)
att1_avg <- mean(att1(D1))

# ============================================================
# Scenario 2: Concave ATT (diminishing marginal harm)
#   D ~ Uniform(5, 20)
#   ATT(d) = 6 * (sqrt(d/5) - 1)   [zero at d=5, concave]
#
#   TWFE overweights high-dose units (large D - D_bar > 0)
#   where the dose-response curve has flattened, and
#   underweights/negatively-weights low-dose units where
#   the slope is steepest. Result: biased causal response.
# ============================================================
D2   <- runif(N, 5, 20)
att2 <- function(d) 6 * (sqrt(d / 5) - 1)

p2       <- make_panel(D2, att2)
beta2    <- twfe_coef(p2)
att2_avg <- mean(att2(D2))

# True ACRT at the median dose (derivative of ATT(d) at d=12.5)
acrt2_true <- function(d) 6 * (1 / (2 * sqrt(d / 5))) * (1/5)
acrt2_at_median <- acrt2_true(median(D2))

# ============================================================
# Scenario 3: WRONG SIGN
#   D: 50% at ~5 ug/m3 (rural, low-PM2.5 areas)
#      50% at ~17 ug/m3 (industrial, high-PM2.5 areas)
#   ATT(d) = 20 * exp(-d/4)   [large harm at low dose,
#                               near-zero harm at high dose]
#
#   Motivation: industrial/urban areas may have more
#   pollution-tolerant (selected) populations, or confounders
#   that attenuate the health signal. The DOSE-RESPONSE IS
#   POSITIVE EVERYWHERE — pollution always harms health.
#
#   True avg ATT ≈ 3.0  (positive harm)
#   TWFE estimate ≈ -0.45  (WRONG SIGN!)
#
#   Why: TWFE weights = (D_i - D_bar).
#     Rural units (D~5):      D - D_bar < 0  -> NEGATIVE weight
#     Industrial units (D~17): D - D_bar > 0  -> POSITIVE weight
#   But rural units have large ATT, industrial units have small ATT.
#   TWFE slaps positive weight on small effects, negative weight on
#   large effects -> reports a negative coefficient.
# ============================================================
D3   <- c(rnorm(N / 2, mean = 5, sd = 1), rnorm(N / 2, mean = 17, sd = 1))
att3 <- function(d) 20 * exp(-d / 4)

p3       <- make_panel(D3, att3)
beta3    <- twfe_coef(p3)
att3_avg <- mean(att3(D3))

# ============================================================
# Console summary
# ============================================================
cat("======================================================\n")
cat("TWFE with Continuous Treatment: Simulation Summary\n")
cat("======================================================\n\n")

cat("Scenario 1: Linear ATT(d) = 0.5 * d\n")
cat(sprintf("  TWFE estimate:        %6.3f  (slope of DeltaY on D)\n", beta1))
cat(sprintf("  True avg ATT:         %6.3f  (E[ATT(D_i)])\n", att1_avg))
cat(sprintf("  Ratio (true/TWFE):    %6.1fx\n", att1_avg / beta1))
cat("  Problem: different estimands — TWFE gives marginal\n")
cat("           causal response, not average level effect.\n\n")

cat("Scenario 2: Concave ATT(d) = 6*(sqrt(d/5) - 1)\n")
cat(sprintf("  TWFE estimate:        %6.3f\n", beta2))
cat(sprintf("  True avg ATT:         %6.3f\n", att2_avg))
cat(sprintf("  True ACRT at median:  %6.3f\n", acrt2_at_median))
cat("  Problem: variance-weighting overweights high-dose\n")
cat("           units where dose-response is flat.\n\n")

cat("Scenario 3: ATT(d) = 20*exp(-d/4), bimodal D\n")
cat(sprintf("  TWFE estimate:        %6.3f  <- WRONG SIGN!\n", beta3))
cat(sprintf("  True avg ATT:         %6.3f  (positive harm)\n", att3_avg))
cat("  Problem: TWFE negatively weights rural low-dose\n")
cat("           units (large ATT) and positively weights\n")
cat("           industrial high-dose units (small ATT).\n")
cat("           Says PM2.5 is beneficial. It is not.\n\n")

# ============================================================
# Figures
# ============================================================
theme_set(theme_classic(base_size = 11) +
  theme(plot.background = element_blank()))

d_seq <- seq(1, 25, by = 0.1)

# --- Helper: extract first differences from panel ---
get_deltas <- function(panel) {
  panel %>%
    select(unit, period, D, y) %>%
    pivot_wider(names_from = period, values_from = y,
                names_prefix = "y") %>%
    mutate(delta_y = y1 - y0)
}

# ============================================================
# Figure 1: Overview — all three scenarios
# ============================================================

# Dose-response curves
dr_data <- tibble(
  d = d_seq,
  `1: Linear`  = att1(d_seq),
  `2: Concave` = att2(d_seq),
  `3: Declining` = att3(d_seq)
) %>%
  pivot_longer(-d, names_to = "scenario", values_to = "att")

fig1_dr <- dr_data %>%
  ggplot(aes(x = d, y = att, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_brewer(palette = "Set1", name = "") +
  labs(x = "Dose D (µg/m³ PM₂.₅)", y = "ATT(d|d)",
       title = "True dose-response functions") +
  theme(legend.position = "bottom")

# Summary comparison: TWFE vs avg ATT
summary_df <- tibble(
  scenario = rep(c("1: Linear", "2: Concave", "3: Declining"), 2),
  estimator = rep(c("TWFE estimate", "True avg ATT"), each = 3),
  value = c(beta1, beta2, beta3, att1_avg, att2_avg, att3_avg)
)

fig1_bars <- summary_df %>%
  ggplot(aes(x = scenario, y = value, fill = estimator)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(
    values = c("TWFE estimate" = "#E41A1C", "True avg ATT" = "#377EB8"),
    name = ""
  ) +
  labs(x = NULL, y = "Estimate",
       title = "TWFE vs truth across scenarios") +
  theme(legend.position = "bottom")

fig1 <- fig1_dr | fig1_bars
ggsave("figures/fig-twfe-continuous-overview.png", fig1,
       width = 10, height = 4.5, dpi = 150)

# ============================================================
# Figure 2: Deep dive on Scenario 3 (wrong sign)
# ============================================================
deltas3 <- get_deltas(p3)
D3_bar  <- mean(D3)

# Panel A: Bimodal distribution of D
pA <- tibble(D = D3) %>%
  ggplot(aes(x = D)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white",
                 linewidth = 0.2) +
  geom_vline(xintercept = D3_bar, linetype = "dashed", color = "red",
             linewidth = 0.8) +
  annotate("text", x = D3_bar + 0.6, y = 75,
           label = sprintf("D̄ = %.1f", D3_bar),
           color = "red", hjust = 0, size = 3.5) +
  labs(x = "PM₂.₅ (µg/m³)", y = "Count",
       title = "A. Bimodal exposure distribution",
       subtitle = "50% rural (D ≈ 5), 50% industrial (D ≈ 17)")

# Panel B: True ATT(d) — positive everywhere
pB <- tibble(d = d_seq, att = att3(d_seq)) %>%
  ggplot(aes(x = d, y = att)) +
  geom_line(color = "#377EB8", linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = att3_avg, linetype = "dotted",
             color = "#377EB8", linewidth = 0.8) +
  annotate("text", x = 20, y = att3_avg + 0.8,
           label = sprintf("E[ATT] = %.2f", att3_avg),
           color = "#377EB8", size = 3.5) +
  labs(x = "PM₂.₅ dose (µg/m³)", y = "ATT(d|d)",
       title = "B. True dose-response",
       subtitle = "Positive everywhere — pollution always harmful")

# Panel C: TWFE weights (D_i - D_bar)
pC <- deltas3 %>%
  mutate(weight = D - D3_bar,
         sign   = if_else(weight < 0, "Negative (below D̄)",
                                      "Positive (above D̄)")) %>%
  ggplot(aes(x = D, y = weight, color = sign)) +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_vline(xintercept = D3_bar, linetype = "dashed",
             color = "grey40", linewidth = 0.6) +
  scale_color_manual(
    values = c("Negative (below D̄)" = "#E41A1C",
               "Positive (above D̄)" = "#377EB8"),
    name = "TWFE weight"
  ) +
  labs(x = "PM₂.₅ (µg/m³)", y = "D_i − D̄",
       title = "C. TWFE weights by unit",
       subtitle = "Rural low-dose units (large ATT) get negative weight") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9))

# Panel D: Scatter ΔY vs D with TWFE line
pD <- deltas3 %>%
  ggplot(aes(x = D, y = delta_y)) +
  geom_point(alpha = 0.2, size = 0.6, color = "grey60") +
  geom_smooth(method = "lm", color = "#E41A1C", se = FALSE,
              linewidth = 1.1) +
  geom_hline(yintercept = att3_avg, linetype = "dotted",
             color = "#377EB8", linewidth = 0.8) +
  annotate("text", x = 21, y = att3_avg + 0.8,
           label = sprintf("True E[ATT] = %.2f", att3_avg),
           color = "#377EB8", size = 3.2, hjust = 1) +
  annotate("text", x = 3, y = min(deltas3$delta_y) + 2,
           label = sprintf("TWFE slope = %.2f", beta3),
           color = "#E41A1C", size = 3.2, hjust = 0) +
  labs(x = "PM₂.₅ (µg/m³)", y = "ΔY (post − pre)",
       title = "D. TWFE regression of ΔY on D",
       subtitle = "Negative slope despite positive effects everywhere")

fig2 <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Scenario 3: TWFE Gets the Wrong Sign",
    subtitle = paste0(
      "True avg ATT = +", round(att3_avg, 2),
      " (PM₂.₅ is harmful),  TWFE estimate = ",
      round(beta3, 2), " (says it's beneficial)"
    )
  )

ggsave("figures/fig-twfe-continuous-wrongsign.png", fig2,
       width = 10, height = 8, dpi = 150)

# ============================================================
# Figure 3: Weight anatomy for Scenario 3
#   Shows each unit's ATT(D_i) vs its TWFE weight (D_i - D_bar)
#   Negative correlation = the units with large effects get
#   negatively weighted => net estimate is negative
# ============================================================
weight_df <- tibble(
  D       = D3,
  att     = att3(D3),
  weight  = D3 - D3_bar,
  sign    = if_else(D3 < D3_bar, "Negative weight", "Positive weight")
)

fig3 <- weight_df %>%
  ggplot(aes(x = weight, y = att, color = sign)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c("Negative weight" = "#E41A1C", "Positive weight" = "#377EB8"),
    name = ""
  ) +
  annotate("text", x = -10, y = 12,
           label = "Large ATT,\nnegative weight\n(rural areas)",
           color = "#E41A1C", size = 3.5) +
  annotate("text", x = 5, y = 0.3,
           label = "Small ATT,\npositive weight\n(industrial areas)",
           color = "#377EB8", size = 3.5) +
  labs(
    x = "TWFE weight (D_i − D̄)",
    y = "True ATT(D_i)",
    title = "TWFE weight vs. unit's true treatment effect",
    subtitle = paste0(
      "Correlation = ", round(cor(weight_df$weight, weight_df$att), 3),
      "  —  units with largest effects get negatively weighted"
    )
  ) +
  theme(legend.position = "none")

ggsave("figures/fig-twfe-weight-anatomy.png", fig3,
       width = 7, height = 5, dpi = 150)

cat("Figures saved to figures/\n")

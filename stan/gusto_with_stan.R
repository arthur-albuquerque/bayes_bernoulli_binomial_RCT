# Ensures the package "pacman" is installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(Hmisc, tidyverse, brms, tidybayes)

# Load data
load(url('http://hbiostat.org/data/gusto.rda'))


# Data preparation
gusto <- 
  upData(gusto, keep=Cs(day30, tx, age, Killip, sysbp, pulse, pmi, miloc)) %>%
  as_tibble() %>%
  mutate(
    age_s = (age - mean(age)) / (sd(age)),
    sbp_s = (sysbp - mean(sysbp)) / (sd(sysbp)),
    pulse_s = (pulse - mean(pulse)) / (sd(pulse)),
    id = seq(1:nrow(gusto)),
    id = factor(id),
    age_c = cut(age_s, breaks = c(-Inf, 0, Inf), labels = c("low", "high")),
    pulse_c = cut(pulse_s, breaks = c(-Inf, 0, Inf), labels = c("low", "high")),
    sbp_c = cut(sbp_s, breaks = c(-Inf, 0, Inf), labels = c("low", "high")),
    killip_c = case_when(Killip %in% c("I", "II") ~ "low",
                         Killip %in% c("III", "IV") ~ "high")
  ) %>%
  mutate(
    group_id = group_indices(., age_c, killip_c, pmi),
    group_id = factor(group_id)
  )

# Fit model
m_cont <-
  brm(
    data = gusto,
    family = bernoulli,
    formula =
      day30 ~ 0 + Intercept + tx + age_s + pulse_s + sbp_s + 
      Killip + pmi + miloc,
    prior = c(prior(normal(-2.5, 0.75), class = b, coef = "Intercept"),
              prior(normal(0, 0.5), class = b, coef = "txSK"),
              prior(normal(0, 0.5), class = b, coef = "txtPA"),
              prior(normal(0, 0.5), class = b, coef = "age_s"),
              prior(normal(0, 0.5), class = b, coef = "pulse_s"),
              prior(normal(0, 0.5), class = b, coef = "sbp_s"),
              prior(normal(0, 0.5), class = b, coef = "KillipII"),
              prior(normal(0, 0.5), class = b, coef = "KillipIII"),
              prior(normal(0, 0.5), class = b, coef = "KillipIV"),
              prior(normal(0, 0.5), class = b, coef = "milocAnterior"),
              prior(normal(0, 0.5), class = b, coef = "milocOther"),
              prior(normal(0, 0.5), class = b, coef = "pmiyes")),
    iter = 2000, warmup = 1000, 
    chains = 4, cores = 4,
    seed = 11,
    backend = "cmdstanr"
  )

m <- cmdstanr::cmdstan_model("gusto.stan")
X <- model.matrix(
  day30 ~ 1 + tx + age_s + pulse_s + sbp_s + 
    Killip + pmi + miloc, data = gusto
) %>% as.data.frame() %>% 
  rename(Intercept := `(Intercept)`) %>% 
  as.matrix()
data <- list(
  N = nrow(gusto),
  Y = gusto$day30,
  X  = X,
  K = 12,
  prior_only = 0,
  n_risks = 50
)
fit <- m$sample(data = data, seed = 123, parallel_chains = 4)
draws <- spread_draws(fit, baseline_events[i], treated_events[i]) %>% 
  mutate(delta_events = baseline_events - treated_events)

p_avoid <- draws %>% 
  group_by(i) %>% 
  summarise(d = mean(delta_events >= 1))

plot_by_risk <- function(risk) {
  draws %>% 
    pivot_longer(cols = contains('events')) %>% 
    mutate(
      name = factor(
        name %>% str_replace("_", " "),
        paste0(
          c("baseline", "treated", "delta"), 
          " events"
        )
      )
    ) %>% 
    filter(i == risk) %>% 
    ggplot(aes(x = value, y = name)) +
    stat_halfeye() +
    geom_vline(xintercept = 0, linetype = 2) +
    theme_ggdist() +
    scale_x_continuous(breaks = scales::pretty_breaks(10)) +
    labs(x = "Events per 1000 patients", y = NULL,
         title = str_glue("Baseline risk: {risk}%"),
         subtitle = str_glue("P(Avoid at least one event) = {p_avoid[risk, 'd']*100}%"),
         caption = "Delta is baseline - treated")
}

plot_by_risk(1)
plot_by_risk(10)
plot_by_risk(20)
plot_by_risk(30)

p_avoid %>%  
  ggplot(aes(i, d)) +
  geom_line() +
  scale_y_continuous(labels = scales::percent,
                     breaks = scales::pretty_breaks(10)) +
  labs(x = "Baseline risk", 
       subtitle = "Probability of avoiding at least one event every 1000 treated patients",
       title = "P(Baseline events - Treated events \u2265 1)",
       y = "Probability") +
  theme_bw()

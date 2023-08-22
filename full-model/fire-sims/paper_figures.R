library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)

# following code is for extracting from the actual model fit -------
burn_fits <- paste0("full-model/fire-sims/joint/sigma-ri/csv-fits/",
                    list.files(path = "full-model/fire-sims/joint/sigma-ri/csv-fits",
                               pattern = ".csv", recursive = TRUE))
best_fit <- burn_fits[grepl("theta-ri_gamma-ri", burn_fits)]
nfits <- length(burn_fits)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- burn_fits[(3*i-2):(3*i)]
}

best_group <- list(best_fit)
burn_names <- lapply(fit_groups, function(x) str_remove(str_remove(basename(x[1]), "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv"), "joint_")) %>% unlist()
# 
extraction <- function(file_group, burn_name) {
  object <- as_cmdstan_fit(file_group)
  # betas <- object$draws(variables = "beta")
  ri_matrix <- object$draws(variables = "ri_matrix")
  reg <- object$draws(variables = "reg")
  temp <- list(reg, ri_matrix)
  names(temp) <- c("reg", "ri_matrix")
  assign(burn_name, temp, parent.frame())
  rm(object)
  gc()
}

extraction(best_group[[1]], "sigma_ri_theta-ri_gamma-ri")

pegpd <- function(y, kappa, sigma, xi) (1 - (1 + xi * (y/sigma))^(-1/xi))^kappa
qegpd <- function(p, kappa, sigma, xi) (sigma/xi) * ( (1 - p^(1/kappa) )^-xi - 1)
rlevel <- function(m, kappa, sigma, xi) {
  p <- (m-1)/m * (1-pegpd(1.001, kappa, sigma, xi)) + pegpd(1.001, kappa, sigma, xi)
  return(qegpd(p, kappa, sigma, xi))
}

kappa_vals <- `sigma_ri_theta-ri_gamma-ri`$reg %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>% 
  separate_wider_delim(cols = "name", delim = ",", names = c("time", "region")) %>%
  mutate(time = as.numeric(gsub("reg\\[", "", time)),
         region = as.numeric(gsub("\\]", "", region)),
         kappa = exp(value)) %>% select(-value)
  
rand_int <- `sigma_ri_theta-ri_gamma-ri`$ri_matrix %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>% 
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "time", "region"))

rand_int <- rand_int %>% select(-time) %>% distinct()
rand_int <- rand_int %>% 
  mutate(param = as.numeric(gsub("ri_matrix\\[", "", param)),
         region = as.numeric(gsub("\\]", "", region)))
rand_int <- rand_int %>%
  mutate(param = as.character(param),
         param = case_when(param == "1" ~ "sigma",
                           param == "2" ~ "xi",
                           TRUE ~ param),
         value = exp(value)) %>%
  pivot_wider(names_from = param, values_from = value)

all_params <- kappa_vals %>% left_join(rand_int)
returns <- all_params %>% mutate(yr50 = rlevel(50, kappa, sigma, xi))
returns_summary <- returns %>% group_by(time, region) %>% 
  summarize(mean50 = mean(yr50[is.finite(yr50)]), med50 = median(yr50[is.finite(yr50)]), 
            lower25 = quantile(yr50[is.finite(yr50)], probs = 0.025), upper975 = quantile(yr50[is.finite(yr50)], probs = 0.975),
            lower50 = quantile(yr50[is.finite(yr50)], probs = 0.05), upper95 = quantile(yr50[is.finite(yr50)], probs = 0.95)) %>%
  ungroup()

region_key <- readRDS(file = "./full-model/data/processed/region_key.rds")
full_reg_key <- as_tibble(region_key) %>% 
  mutate(region = c(1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE))
reg_cols <- full_reg_key$region

returns_regional <- returns_summary %>% left_join(full_reg_key)
 
p <- returns_regional %>% ggplot() + 
  geom_line(aes(x=time, y=med50, group = region, alpha = 0.5), linewidth = 0.2, color = 'grey') + 
  geom_ribbon(aes(x=time, ymin=lower25, ymax=upper975, group = region, fill = NA_L1CODE, alpha = 0.05)) +
  facet_wrap(. ~ NA_L1CODE, nrow = 2) + theme_classic() + theme(legend.position = "none")

file_name <- "full-model/figures/paper/50yr_returns_bylevel1_95ci.png"
ggsave(file_name, p, dpi = 320, bg = "white")

p <- returns_regional %>% ggplot() + 
  geom_line(aes(x=time, y=med50, group = region, alpha = 0.5), linewidth = 0.2, color = 'grey') + 
  geom_ribbon(aes(x=time, ymin=lower50, ymax=upper95, group = region, fill = NA_L1CODE, alpha = 0.05)) +
  facet_wrap(. ~ NA_L1CODE, nrow = 2) + theme_classic() + theme(legend.position = "none")

file_name <- "full-model/figures/paper/50yr_returns_bylevel1_90ci.png"
ggsave(file_name, p, dpi = 320, bg = "white")

ri_map <- `sigma_ri_theta-ri_gamma-ri`$ri_matrix %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>% 
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "time", "region")) %>% select(-time) %>% distinct() %>%
  mutate(param = as.character(gsub("ri_matrix\\[", "", param)),
         region = as.numeric(gsub("\\]", "", region)),
         param = case_when(param == "1" ~ "sigma",
                           param == "2" ~ "xi",
                           TRUE ~ param))

ri_map <- ri_map %>% group_by(param, region) %>% summarize(med_val = median(value))
ri_map <- ri_map %>% ungroup()
xi_map <- ri_map %>% filter(param == "xi") %>% select(-param) %>% mutate(med_val = exp(med_val))
sigma_map <- ri_map %>% filter(param == "sigma") %>% select(-param) %>% mutate(med_val = exp(med_val))

xi_vals_regions <- xi_map %>% left_join(full_reg_key)
eco_xi <- ecoregions_geom %>% 
  mutate(NA_L2CODE = as.factor(NA_L2CODE), 
         NA_L1CODE = as.factor(NA_L1CODE), 
         NA_L3CODE = as.factor(NA_L3CODE)) %>% 
  left_join(xi_vals_regions)

breaks <- classIntervals(c(min(eco_xi$med_val) - .00001, eco_xi$med_val), style = 'fixed', 
                         fixedBreaks = c(0, 0.2, 0.4, 0.6, 0.8, 2.0), intervalClosure = 'left')

eco_xi_cat <- eco_xi %>% 
  mutate(xi_cat = cut(med_val, unique(breaks$brks)))

p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_xi_cat,
          aes(fill=xi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_void() + scale_fill_brewer(palette = 'YlOrRd')
ggsave("full-model/figures/paper/xi_map.png", dpi = 320, bg ='white')


sigma_vals_regions <- sigma_map %>% left_join(full_reg_key)
eco_sigma <- ecoregions_geom %>% 
  mutate(NA_L2CODE = as.factor(NA_L2CODE), 
         NA_L1CODE = as.factor(NA_L1CODE), 
         NA_L3CODE = as.factor(NA_L3CODE)) %>% 
  left_join(sigma_vals_regions)

breaks <- classIntervals(c(min(eco_sigma$med_val) - .00001, eco_sigma$med_val), style = 'fisher', n=5, intervalClosure = 'left')
eco_sigma_cat <- eco_sigma %>% 
  mutate(sigma_cat = cut(med_val, unique(breaks$brks)))

p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_sigma_cat,
          aes(fill=sigma_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_void() + scale_fill_brewer(palette = 'YlOrRd')
ggsave("full-model/figures/paper/sigma_map.png", dpi = 320, bg ='white')



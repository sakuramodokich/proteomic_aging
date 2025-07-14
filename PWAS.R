library(tidyverse)
library(furrr)
library(lme4)
options(future.globals.maxSize=12*1024^3)

dat <- data.table::fread("main_aging_dataset.csv") %>% left_join(data.table::fread("main_aging_cov.csv"),"eid") %>% as_tibble
pro <- data.table::fread("/mnt/vol2/UKB_dataset/asset/olink_protein_wide2923_QC_knn_impute_INT.tsv") %>% as_tibble

dat <- dat %>% semi_join(pro, "eid")
pro <- pro %>% semi_join(dat, "eid")

prot_list <- names(pro)[-1]
aging_list <- c('kdm_delta','phenoage_delta','ltl','frailty_score','healthspan')
cov_list <- c('age','sex','tdi','ethnic','qualification','diet_score','smoke_status','alcohol_freq','bmi')

dat <- dat %>% mutate(
  qualification = factor(qualification),
  smoke_status = factor(smoke_status),
  alcohol_freq = factor(alcohol_freq),
  ass_center = factor(ass_center)
)

multi.regress.tidy <- function(data, omic, exposure, outcome, covs){
  
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  
  pheno <- omic %>% select(eid, all_of(exposure)) %>% dplyr::rename("x" = exposure)
  data <- data %>% dplyr::rename("y"=outcome) %>% left_join(pheno, "eid") %>% drop_na(x, y, covs)
  data$y <- RNOmni::RankNorm(data$y)
  formula <- paste0("y ~ x + ", paste0(covs, collapse = " + "))
  
  fit <- lm(as.formula(formula), data = data)
  res <- fit %>% broom::tidy(conf.int=T) %>% .[2, ]
  res$term <- exposure
  res$outcome <- outcome
  res$N <- nrow(data)
  return(res)
}

multi.regress.mixed.tidy <- function(data, omic, exposure, outcome, covs, random_effects){
  
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  
  pheno <- omic %>% select(eid, all_of(exposure)) %>% dplyr::rename("x" = exposure)
  data <- data %>% dplyr::rename("y"=outcome) %>% left_join(pheno, "eid") %>% drop_na(x, y, covs)
  data$y <- RNOmni::RankNorm(data$y)
  
  fixed_formula <- paste0("y ~ x + ", paste0(covs, collapse = " + "))
  formula <- paste0(fixed_formula, " + ", random_effects)
  
  fit <- lme4::lmer(as.formula(formula), data = data)
  res <- broom.mixed::tidy(fit, conf.int=T, effects = "fixed") %>% filter(term == "x")
  
  res$term <- exposure
  res$outcome <- outcome
  res$N <- nrow(data)
  return(res)
}

full_prot_res <- NULL
for(i in 1:length(aging_list)) {
  plan(multisession, workers = 90)
  print(paste0("Running: ",aging_list[i]))
  prot_res <- future_map_dfr(prot_list, ~ multi.regress.mixed.tidy(dat, pro, ., aging_list[i], cov_list, "(1|ass_center)"), .options = furrr_options(seed = TRUE))
  plan(sequential)
  full_prot_res <- rbind(full_prot_res, prot_res)
}
full_prot_res %>% write_csv("result/step1_PWAS_lmer.csv")

full_prot_res <- NULL
for(i in 1:length(aging_list)) {
  plan(multisession, workers = 90)
  print(paste0("Running: ",aging_list[i]))
  prot_res <- future_map_dfr(prot_list, ~ multi.regress.tidy(dat, pro, ., aging_list[i], cov_list), .options = furrr_options(seed = TRUE))
  prot_res$P_BON <- p.adjust(prot_res$p.value, method="bonferroni")
  plan(sequential)
  full_prot_res <- rbind(full_prot_res, prot_res)
}

full_prot_res %>% write_csv("result/step1_PWAS.csv")

plan(multisession, workers = 90)
prot_res <- future_map_dfr(prot_list, ~ multi.regress.tidy(dat, pro, ., "age", cov_list[cov_list!="age"]), .options = furrr_options(seed = TRUE))
prot_res$P_BON <- p.adjust(prot_res$p.value, method="bonferroni")
plan(sequential)
prot_res %>% write_csv("result/step1_PWAS_age.csv")

# Sex
multi.regress.int.tidy <- function(data, omic, exposure, outcome, covs, interaction_term = NULL){
  
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  
  pheno <- omic %>% select(eid, all_of(exposure)) %>% dplyr::rename("x" = exposure)
  data <- data %>% dplyr::rename("y"=outcome) %>% left_join(pheno, "eid") %>% drop_na(x, y, covs)
  data$y <- RNOmni::RankNorm(data$y)
  formula <- paste0("y ~ x + ", paste0(covs, collapse = " + "))
  
  data <- data %>% dplyr::rename("int" = interaction_term)
  formula <- paste0("y ~ x * int + ", paste0(covs[covs!=interaction_term], collapse = " + "))
  fit <- lm(as.formula(formula), data = data)
  res <- fit %>% broom::tidy(conf.int=T) %>% tail(1)
  
  res$term <- exposure
  res$outcome <- outcome
  res$N <- nrow(data)
  return(res)
}

full_prot_res <- NULL
for(i in 1:length(aging_list)) {
  plan(multisession, workers = 100)
  print(paste0("Running: ",aging_list[i]))
  prot_res <- future_map_dfr(prot_list, ~ multi.regress.tidy(dat %>% filter(sex == 0), pro, ., aging_list[i], cov_list[cov_list!="sex"]), .options = furrr_options(seed = TRUE))
  prot_res$P_FDR <- p.adjust(prot_res$p.value, method="BH")
  prot_res$P_BON <- p.adjust(prot_res$p.value, method="bonferroni")
  plan(sequential)
  full_prot_res <- rbind(full_prot_res, prot_res)
}
full_prot_res$sex_strata <- "female"
full_prot_res2 <- NULL
for(i in 1:length(aging_list)) {
  plan(multisession, workers = 100)
  print(paste0("Running: ",aging_list[i]))
  prot_res <- future_map_dfr(prot_list, ~ multi.regress.tidy(dat %>% filter(sex == 1), pro, ., aging_list[i], cov_list[cov_list!="sex"]), .options = furrr_options(seed = TRUE))
  prot_res$P_FDR <- p.adjust(prot_res$p.value, method="BH")
  prot_res$P_BON <- p.adjust(prot_res$p.value, method="bonferroni")
  plan(sequential)
  full_prot_res2 <- rbind(full_prot_res2, prot_res)
}
full_prot_res2$sex_strata <- "male"
full_prot_res <- rbind(full_prot_res, full_prot_res2)
full_prot_res %>% write_csv("result/PWAS_supp_sex_strata.csv")

# Sex interaction
full_prot_res <- NULL
for(i in 1:length(aging_list)) {
  plan(multisession, workers = 100)
  print(paste0("Running: ",aging_list[i]))
  prot_res <- future_map_dfr(prot_list, ~ multi.regress.int.tidy(dat, pro, ., aging_list[i], cov_list, interaction_term = "sex"), .options = furrr_options(seed = TRUE))
  prot_res$P_FDR <- p.adjust(prot_res$p.value, method="BH")
  prot_res$P_BON <- p.adjust(prot_res$p.value, method="bonferroni")
  plan(sequential)
  full_prot_res <- rbind(full_prot_res, prot_res)
}
full_prot_res %>% write_csv("result/PWAS_supp_sex_interaction.csv")

# Age
multi.regress.int.tidy <- function(data, omic, exposure, outcome, covs, interaction_term = NULL){
  
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  
  pheno <- omic %>% select(eid, all_of(exposure)) %>% dplyr::rename("x" = exposure)
  data <- data %>% dplyr::rename("y"=outcome) %>% left_join(pheno, "eid") %>% drop_na(x, y, covs)
  data$y <- RNOmni::RankNorm(data$y)
  formula <- paste0("y ~ x + ", paste0(covs, collapse = " + "))
  
  data <- data %>% dplyr::rename("int" = interaction_term)
  formula <- paste0("y ~ x * int + ", paste0(covs[covs!=interaction_term], collapse = " + "))
  fit <- lm(as.formula(formula), data = data)
  res <- fit %>% broom::tidy(conf.int=T) %>% tail(1)
  
  res$term <- exposure
  res$outcome <- outcome
  res$N <- nrow(data)
  return(res)
}
dat$age <- ifelse(dat$age >= 65, 1, 0) # 1:>=65

full_prot_res <- NULL
for(i in 1:length(aging_list)) {
  plan(multisession, workers = 100)
  print(paste0("Running: ",aging_list[i]))
  prot_res <- future_map_dfr(prot_list, ~ multi.regress.tidy(dat %>% filter(age == 0), pro, ., aging_list[i], cov_list[cov_list!="age"]), .options = furrr_options(seed = TRUE))
  prot_res$P_FDR <- p.adjust(prot_res$p.value, method="BH")
  prot_res$P_BON <- p.adjust(prot_res$p.value, method="bonferroni")
  plan(sequential)
  full_prot_res <- rbind(full_prot_res, prot_res)
}
full_prot_res$sex_strata <- "<65"
full_prot_res2 <- NULL
for(i in 1:length(aging_list)) {
  plan(multisession, workers = 100)
  print(paste0("Running: ",aging_list[i]))
  prot_res <- future_map_dfr(prot_list, ~ multi.regress.tidy(dat %>% filter(age == 1), pro, ., aging_list[i], cov_list[cov_list!="age"]), .options = furrr_options(seed = TRUE))
  prot_res$P_FDR <- p.adjust(prot_res$p.value, method="BH")
  prot_res$P_BON <- p.adjust(prot_res$p.value, method="bonferroni")
  plan(sequential)
  full_prot_res2 <- rbind(full_prot_res2, prot_res)
}
full_prot_res2$sex_strata <- ">=65"
full_prot_res <- rbind(full_prot_res, full_prot_res2)
full_prot_res %>% write_csv("result/PWAS_supp_age_strata.csv")

# Age interaction
full_prot_res <- NULL
for(i in 1:length(aging_list)) {
  plan(multisession, workers = 100)
  print(paste0("Running: ",aging_list[i]))
  prot_res <- future_map_dfr(prot_list, ~ multi.regress.int.tidy(dat, pro, ., aging_list[i], cov_list, interaction_term = "age"), .options = furrr_options(seed = TRUE))
  prot_res$P_FDR <- p.adjust(prot_res$p.value, method="BH")
  prot_res$P_BON <- p.adjust(prot_res$p.value, method="bonferroni")
  plan(sequential)
  full_prot_res <- rbind(full_prot_res, prot_res)
}
full_prot_res %>% write_csv("result/PWAS_supp_age_interaction.csv")

# Age distribution
library(ggridges)
library(ggstatsplot)
library(ggside)
library(patchwork)

p1 <- ggscatterstats(data = dat, x = age, y = kdm_delta, type = "n", marginal=F, results.subtitle = F,
                     point.args = list(alpha = 0.01, size = 3, color = '#2078b4'), 
                     smooth.line.args = list(method = "lm", color="#e35f2b", fill = "#e35f2b", alpha=0.15, size=0.9)) +
  geom_xsidedensity(aes(y = after_stat(density)), fill='#7aa4cc') +
  geom_ysidedensity(aes(x = after_stat(density)), fill='#7aa4cc') +
  ggpubr::stat_cor(method = "spearman", size = 4, cor.coef.name = "r") +
  xlab("Age") + ylab("KDM-BA acceleration") + theme_classic() +
  theme(axis.title = element_text(size=12, color = "black"),
        axis.text = element_text(size=10, color = "black"),
        ggside.panel.scale = 0.25, 
        axis.text.x = element_text(angle = 90, vjust = 0.5))
p2 <- ggscatterstats(data = dat, x = age, y = phenoage_delta, type = "n", marginal=F, results.subtitle = F,
                     point.args = list(alpha = 0.01, size = 3, color = '#2078b4'), 
                     smooth.line.args = list(method = "lm", color="#e35f2b", fill = "#e35f2b", alpha=0.15, size=0.9)) +
  geom_xsidedensity(aes(y = after_stat(density)), fill='#7aa4cc') +
  geom_ysidedensity(aes(x = after_stat(density)), fill='#7aa4cc') +
  ggpubr::stat_cor(method = "spearman", size = 4, cor.coef.name = "r") +
  xlab("Age") + ylab("PhenoAge acceleration") + theme_classic() +
  theme(axis.title = element_text(size=12, color = "black"),
        axis.text = element_text(size=10, color = "black"),
        ggside.panel.scale = 0.25,
        axis.text.x = element_text(angle = 90, vjust = 0.5))
p3 <- ggscatterstats(data = dat, x = age, y = frailty_score, type = "n", marginal=F, results.subtitle = F,
                     point.args = list(alpha = 0.01, size = 3, color = '#2078b4'), 
                     smooth.line.args = list(method = "lm", color="#e35f2b", fill = "#e35f2b", alpha=0.15, size=0.9)) +
  geom_xsidedensity(aes(y = after_stat(density)), fill='#7aa4cc') +
  geom_ysidedensity(aes(x = after_stat(density)), fill='#7aa4cc') +
  ggpubr::stat_cor(method = "spearman", size = 4, cor.coef.name = "r") +
  xlab("Age") + ylab("Frailty index") + theme_classic() +
  theme(axis.title = element_text(size=12, color = "black"),
        axis.text = element_text(size=10, color = "black"),
        ggside.panel.scale = 0.25,
        axis.text.x = element_text(angle = 90, vjust = 0.5))
p4 <- ggscatterstats(data = dat, x = age, y = healthspan, type = "n", marginal=F, results.subtitle = F,
                     point.args = list(alpha = 0.01, size = 3, color = '#2078b4'), 
                     smooth.line.args = list(method = "lm", color="#e35f2b", fill = "#e35f2b", alpha=0.15, size=0.9)) +
  geom_xsidedensity(aes(y = after_stat(density)), fill='#7aa4cc') +
  geom_ysidedensity(aes(x = after_stat(density)), fill='#7aa4cc') +
  ggpubr::stat_cor(method = "spearman", size = 4, cor.coef.name = "r") +
  xlab("Age") + ylab("Healthspan") + theme_classic() +
  theme(axis.title = element_text(size=12, color = "black"),
        axis.text = element_text(size=10, color = "black"),
        ggside.panel.scale = 0.25,
        axis.text.x = element_text(angle = 90, vjust = 0.5))
p5 <- ggscatterstats(data = dat, x = age, y = ltl, type = "n", marginal=F, results.subtitle = F,
                     point.args = list(alpha = 0.01, size = 3, color = '#2078b4'), 
                     smooth.line.args = list(method = "lm", color="#e35f2b", fill = "#e35f2b", alpha=0.15, size=0.9)) +
  geom_xsidedensity(aes(y = after_stat(density)), fill='#7aa4cc') +
  geom_ysidedensity(aes(x = after_stat(density)), fill='#7aa4cc') +
  ggpubr::stat_cor(method = "spearman", size = 4, cor.coef.name = "r") +
  xlab("Age") + ylab("Leukocyte telomere length") + theme_classic() +
  theme(axis.title = element_text(size=12, color = "black"),
        axis.text = element_text(size=10, color = "black"),
        ggside.panel.scale = 0.25,
        axis.text.x = element_text(angle = 90, vjust = 0.5))

data <- dat %>%
  mutate(age_group = cut(age, 
                         breaks = seq(40, 75, by = 5),
                         labels = paste0(seq(40, 70, by = 5), "-", seq(44, 74, by = 5)),
                         right = FALSE)) %>% drop_na(age_group)

long_data <- data %>%
  pivot_longer(
    cols = c("frailty_score", "healthspan", "ltl", "phenoage_delta", "kdm_delta"),
    names_to = "aging_phenotype",
    values_to = "value"
  )
long_data$aging_phenotype <- factor(
  long_data$aging_phenotype,
  levels = c("kdm_delta", "phenoage_delta", "frailty_score", "healthspan", "ltl"),
  labels = c("KDM-BA acceleration", "PhenoAge acceleration", "Frailty index", "Healthspan", "Leukocyte telomere length")
)
p6 <- ggplot(long_data, aes(x = value, y = age_group, fill = age_group)) +
  facet_wrap(~ aging_phenotype, scales = "free_x", nrow=1) +
  geom_density_ridges(alpha = 0.7, scale = 0.9,
                      rel_min_height = 0.01,
                      alpha = 0.7,
                      quantile_lines = TRUE,
                      quantiles = 2 ) +
  scale_fill_viridis_d() +
  labs(x = NULL, y = "Age groups") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size=12, color = "black"),
        axis.title = element_text(size=12, color = "black"),
        axis.text = element_text(size=10, color = "black"))
design <- "
  ABC
  DE#
  FFF
"
p <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(design = design)
ggsave("supp1.png", plot = p, width = 16, height = 13)

plot_cor <- function(trait, data, xt, x, y) {
  data %>% filter(outcome == trait) %>% ggplot(aes(x = zscore1, y = zscore2)) +
    geom_smooth(method = "lm", color = "#8c8c8c", se = F) +
    ggrastr::geom_point_rast(alpha = 0.2, color = "#1f77b4", stroke = 1) +
    ggprism::theme_prism(base_fontface = "plain", base_size = 4) +
    ggpubr::stat_cor(method = "spearman", size = 4, cor.coef.name = "r") +
    labs(y = y, x = x) + ggtitle(xt) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'),
      plot.title = element_text(size=12, color = "black"),
      axis.title = element_text(size=12, color = "black"),
      axis.text.x = element_text(size=10, color = "black"),
      axis.text.y = element_text(size=10, color = "black"),
      axis.ticks.length = unit(rep(0.25,1),'lines')
    )
}

plot1 <- data.table::fread("result/step1_PWAS.csv")
plot1 <- plot1 %>% mutate(zscore = estimate / std.error)
plot2 <- data.table::fread("result/step1_PWAS_age.csv")
plot2 <- plot2 %>% mutate(zscore = estimate / std.error)
plot2 <- rbind(plot2, plot2, plot2, plot2, plot2)

plot <- tibble(term=plot1$term, zscore1=plot1$zscore, zscore2=plot2$zscore, outcome = plot1$outcome)
x <- 'Z-score'; y <- 'Z-score (chronological age)'
p1 <- plot_cor('kdm_delta', plot, "KDM-BA acceleration", x = x, y = y)
p2 <- plot_cor('phenoage_delta', plot, "PhenoAge acceleration", x = x, y = y)
p3 <- plot_cor('frailty_score', plot, "Frailty index", x = x, y = y)
p4 <- plot_cor('healthspan', plot, "Healthspan", x = x, y = y)
p5 <- plot_cor('ltl', plot, "Leukocyte telomere length", x = x, y = y)

design <- "
  ABC
  EF#
"
p <- p1 + p2 + p3 + p4 + p5 + plot_layout(design = design) 
ggsave("PWAS_supp.png", plot = p, width = 12, height = 8)

plot1 <- data.table::fread("result/step1_PWAS.csv")
plot1 <- plot1 %>% mutate(zscore = estimate / std.error)
plot2 <- data.table::fread("result/step1_PWAS_lmer.csv")
plot2 <- plot2 %>% mutate(zscore = estimate / std.error)

plot <- tibble(term=plot1$term, zscore1=plot1$zscore, zscore2=plot2$zscore, outcome = plot1$outcome)
x <- 'Z-score (linear model)'; y <- 'Z-score (linear mixed model)'
p1 <- plot_cor('kdm_delta', plot, "KDM-BA acceleration", x = x, y = y)
p2 <- plot_cor('phenoage_delta', plot, "PhenoAge acceleration", x = x, y = y)
p3 <- plot_cor('frailty_score', plot, "Frailty index", x = x, y = y)
p4 <- plot_cor('healthspan', plot, "Healthspan", x = x, y = y)
p5 <- plot_cor('ltl', plot, "Leukocyte telomere length", x = x, y = y)

design <- "
  ABC
  EF#
"
p <- p1 + p2 + p3 + p4 + p5 + plot_layout(design = design) 
ggsave("PWAS_supp2.png", plot = p, width = 12, height = 8)

library(tidyverse)
library(furrr)
library(survival)
options(future.globals.maxSize=12*1024^4)

nmr_scatter_path <- "NMR"
nmr_scatter_list <- list.files(nmr_scatter_path) %>% file.path(nmr_scatter_path,.)

dat <- data.table::fread("main_aging_dataset.csv") %>% left_join(data.table::fread("main_aging_cov.csv"),"eid") %>% as_tibble
sig_pro <- data.table::fread("result/step1_MR.csv") %>% filter(P_BON < 0.05) %>% pull(exposure) %>% unique
sig_pro <- gsub("[-]", "_", sig_pro)
pro <- data.table::fread("olink_protein_knn_impute_INT.tsv", select = c('eid', sig_pro)) %>% as_tibble
names(pro) <- gsub("[_]", "-", names(pro))
sig_pro <- names(pro)[-1]

dat <- dat %>% semi_join(pro, "eid")
pro <- pro %>% semi_join(dat, "eid")
dat <- dat %>% left_join(pro, "eid")

aging_list <- c('kdm_delta','phenoage_delta','ltl','frailty_score','healthspan')
cov_list <- c('age','sex','tdi','ethnic','qualification','diet_score','smoke_status','alcohol_freq','bmi')

dat <- dat %>% mutate(
  qualification = factor(qualification),
  smoke_status = factor(smoke_status),
  alcohol_freq = factor(alcohol_freq)
)

# Protein - Metabolite (LM)
multi.regress.tidy <- function(data, trait, exposure, covariate){
  
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  
  omic <- data.table::fread(trait, header = T) %>% semi_join(data, "eid")
  trait_name <- names(omic)[2]
  names(omic)[2] <- "y"
  data <- data %>% dplyr::select(all_of(c('eid',exposure,covariate))) %>% dplyr::rename("x"=exposure) %>% left_join(omic, "eid") %>% drop_na(x, y, covariate)
  data$y <- RNOmni::RankNorm(data$y)
  formula <- paste0("y ~ x + ", paste0(covariate,collapse = " + "))
  res <- lm(as.formula(formula), data=data) %>% broom::tidy(conf.int=T) %>% .[2, ]
  res$term <- trait_name
  res$exposure <- exposure
  res$N <- nrow(data)
  return(res)
}
multi.regress.flip.tidy <- function(data, trait, exposure, covariate){
  
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  
  omic <- data.table::fread(exposure, header = T) %>% semi_join(data, "eid")
  exposure_name <- names(omic)[2]
  names(omic)[2] <- "x"
  data <- data %>% dplyr::select(all_of(c('eid',trait,covariate))) %>% dplyr::rename("y"=trait) %>% left_join(omic, "eid") %>% drop_na(x, y, covariate)
  data$x <- RNOmni::RankNorm(data$x)
  data$y <- RNOmni::RankNorm(data$y)
  formula <- paste0("y ~ x + ", paste0(covariate,collapse = " + "))
  res <- lm(as.formula(formula), data=data) %>% broom::tidy(conf.int=T) %>% .[2, ]
  res$term <- exposure_name
  res$outcome <- trait
  res$N <- nrow(data)
  return(res)
}


full_nmr_res <- NULL
for (i in 1:length(sig_pro)){
  plan(multisession, workers = 50)
  nmr_res_lm <- future_map_dfr(nmr_scatter_list, ~ multi.regress.tidy(dat, ., sig_pro[i], cov_list), .options = furrr_options(seed = TRUE))
  nmr_res_lm$P_BON <- p.adjust(nmr_res_lm$p.value, method = "bonferroni")
  full_nmr_res <- rbind(full_nmr_res, nmr_res_lm)
  plan(sequential)
}
full_nmr_res %>% write_csv("step5_pair_protein_metabolite.csv")

# Metabolite - Aging (LM)
full_nmr_res <- NULL
for (i in 1:length(aging_list)){
  plan(multisession, workers = 50)
  nmr_res_lm <- future_map_dfr(nmr_scatter_list, ~ multi.regress.flip.tidy(dat, aging_list[i], ., cov_list), .options = furrr_options(seed = TRUE))
  nmr_res_lm$P_BON <- p.adjust(nmr_res_lm$p.value, method = "bonferroni")
  full_nmr_res <- rbind(full_nmr_res, nmr_res_lm)
  plan(sequential)
}
full_nmr_res %>% write_csv("step5_pair_aging_metabolite.csv")

# Mediation
med.func <- function(data, x, m, y, covariate) {
  RhpcBLASctl::blas_set_num_threads(6)
  RhpcBLASctl::omp_set_num_threads(6)
  format <- function(results){
    tmp <- summary(results)
    acme <- tmp$d.avg
    acme.L <- as.numeric(tmp$d.avg.ci)[1]
    acme.H <- as.numeric(tmp$d.avg.ci)[2]
    acme.p <- tmp$d.avg.p
    ade <- tmp$z.avg
    ade.L <- as.numeric(tmp$z.avg.ci)[1]
    ade.H <- as.numeric(tmp$z.avg.ci)[2]
    ade.p <- tmp$z.avg.p
    prop <- tmp$n.avg
    prop.L <- as.numeric(tmp$n.avg.ci)[1]
    prop.H <- as.numeric(tmp$n.avg.ci)[2]
    prop.p <- tmp$n.avg.p
    return(as_tibble(cbind(acme, acme.L, acme.H, acme.p, ade, ade.L, ade.H, ade.p, prop, prop.L, prop.H, prop.p)))
  }
  order_med <- function(indf){
    formula <- paste0("m ~ x + ", paste0(covariate, collapse = " + "))
    model.m <- lm(as.formula(formula), data = indf)
    formula <- paste0("y ~ x + m + ", paste0(covariate, collapse = " + "))
    model.y <- lm(as.formula(formula), data = indf)
    med.out <- mediation::mediate(model.m, model.y, sims = 500, treat = "x", mediator = "m", covariates = covariate)
    med.out %>% format(.) %>% return
  }
  mediator <- data.table::fread(paste0("NMR/",m,".tsv")) %>% as_tibble
  data <- data %>% left_join(mediator, "eid")
  data <- data %>% dplyr::select(all_of(x),all_of(y),all_of(covariate),all_of(m)) %>% 
    dplyr::rename("x"=x,"y"=y,"m"=m)  %>% drop_na() 
  data$y <- RNOmni::RankNorm(data$y)
  data$m <- RNOmni::RankNorm(data$m)
  data %>% order_med(.) %>% mutate(x = x, y = y, m = m) -> res
  res$n <- nrow(data)
  res %>% return()
}

med_list <- data.table::fread("result/step1_MR.csv") %>% filter(P_BON < 0.05)
full.med.result <- NULL
for (i in 1:nrow(med_list)){
  dat.med <- dat %>% select(all_of(c('eid',med_list$exposure[i],med_list$outcome[i],cov_list)))
  nmr_list <- nmr_scatter_list %>% basename %>% gsub(".tsv","",.)
  plan(multisession, workers = 15)
  med.result <- future_map_dfr(nmr_list, ~ med.func(dat.med, med_list$exposure[i], ., med_list$outcome[i], cov_list), .options = furrr_options(seed = TRUE))
  full.med.result <- rbind(full.med.result, med.result)
  plan(sequential)
}
full.med.result %>% write_csv("step5_mediation_aging_metabolite.csv")

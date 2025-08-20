library(tidyverse)
library(furrr)
library(survival)
options(future.globals.maxSize=12*1024^4)

dat <- data.table::fread("main_aging_dataset.csv") %>% left_join(data.table::fread("main_aging_cov.csv"),"eid") %>% as_tibble
data_path <- "UKB_merged_dataset.tsv"
date <- data.table::fread(data_path, select=c('eid', 'p53_i0')) %>% as_tibble
dat <- dat %>% left_join(date, "eid")
icd_epistart <- data.table::fread("icd_epistart.tsv", header = T) %>% semi_join(dat, "eid") %>% keep(~sum(!is.na(.)) >= 200)
annot <- read_csv("Phecode_map_v1_2_icd10cm_beta.csv", col_types = list("phecode"="character")) %>% filter(!duplicated(phecode))

multicox.tidy <- possibly(function(data, trait, exposure, covariates){
  icd <- data.table::fread(paste0('icd_epistart_',trait,".tsv"), header = T) %>% dplyr::rename("y"=trait) %>% semi_join(data, "eid")
  stopifnot(sum(as.numeric(is.na(icd$y)))!=nrow(icd))
  data <- data %>% dplyr::rename("x"=exposure) %>% semi_join(icd, "eid") %>% add_column(icd[,-1]) %>% filter(y>0 | is.na(y)) %>% tidyr::drop_na(x, covariates)
  data$y_s <- ifelse(is.na(data$y), 0, 1)
  data$y[is.na(data$y) & data$p54_i0 %in% c(11003,11023,11022)] <- as.numeric(difftime(data.table::as.IDate("2022-5-31"), data$p53_i0[is.na(data$y) & data$p54_i0 %in% c(11003,11023,11022)], units = "days"))/365.25
  data$y[is.na(data$y) & data$p54_i0 %in% c(11005,11004)] <- as.numeric(difftime(data.table::as.IDate("2022-8-31"), data$p53_i0[is.na(data$y) & data$p54_i0 %in% c(11005,11004)], units = "days"))/365.25
  data$y[is.na(data$y)] <- as.numeric(difftime(data.table::as.IDate("2022-10-31"), data$p53_i0[is.na(data$y)], units = "days"))/365.25
  
  stopifnot(sum(data$y_s)!=0)
  formula <- paste0("Surv(y, y_s) ~ x + ", paste0(covariates, collapse = " + "))
  res <- coxph(as.formula(formula), data=data) %>% broom::tidy(exp=T, conf.int=T) %>% head(1)
  res$term <- trait
  res$N <- nrow(data)
  res$Ncase <- sum(data$y_s)
  res$Ncontrol <- res$N - res$Ncase
  res$description <- annot$phecode_str[annot$phecode == trait]
  res$category <- annot$exclude_name[annot$phecode == trait]
  res$exposure <- exposure
  return(res)
})

phecode_phewas <- names(icd_epistart)[-1]
rm(icd_epistart);gc()

sig_pro <- data.table::fread("result/step1_MR.csv") %>% filter(P_BON < 0.05) %>% pull(exposure) %>% unique
sig_pro <- gsub("[-]", "_", sig_pro)
pro <- data.table::fread("olink_protein_QC_knn_impute_INT.tsv", select = c('eid', sig_pro)) %>% as_tibble
names(pro) <- gsub("[_]", "-", names(pro))
sig_pro <- names(pro)[-1]
cov_list <- c('age','sex','tdi','ethnic','qualification','diet_score','smoke_status','alcohol_freq','bmi')
dat <- dat %>% left_join(pro, "eid")

dat <- dat %>% mutate(
  qualification = factor(qualification),
  smoke_status = factor(smoke_status),
  alcohol_freq = factor(alcohol_freq)
)

full_results <- NULL
chunk_size <- 800
for (k in 1:length(sig_pro)) {
  protein <- sig_pro[k]
  for (i in seq(0, length(phecode_phewas), by=chunk_size)) {
    plan(multisession, workers = 50)
    start <- i+1
    end <- ifelse(i+chunk_size>length(phecode_phewas),length(phecode_phewas),i+chunk_size)
    print(c(start,end))
    results <- future_map(phecode_phewas[start:end], ~ multicox.tidy(dat, ., protein, cov_list), .options = furrr_options(seed = T)) %>% purrr::reduce(rbind)
    full_results <- rbind(full_results, results)
    plan(sequential)
  }
}
names(full_results) <- c('Phecode','HR','SE','T_value','P','CI_low','CI_high','N','Ncase','Ncontrol','Description','Category','Protein')
full_results %>% write_csv("result/step4_phewas_protein.csv")

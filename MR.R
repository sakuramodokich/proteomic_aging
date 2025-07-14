library(tidyverse)
library(TwoSampleMR)
library(furrr)

pQTL <- data.table::fread("/mnt/vol1/Resources/UKB_PPP_pQTL/UKB_PPP_5e8_hg38_Merged_SNP_n2677_Clumping_1MB_0.001.tsv") %>% as_tibble
pQTL_list <- pQTL$PROTEIN %>% unique 
aging_list <- c('kdm_delta','phenoage_delta','ltl','frailty_score','healthspan')
aging_list <- "healthspan_new"

mr_combined_val <- possibly(function(x){
  x2 <- dat_to_MRInput(x)
  if(nrow(x) == 1) {
    res <- mr(x, method_list = 'mr_wald_ratio')
  } else {
    res <- mr(x, method_list = c('mr_ivw','mr_weighted_median','mr_egger_regression','mr_raps'))
  }
  h_pval <- NULL
  if (nrow(x) != 1) {
    try(
      {
        res_presso <- MRPRESSO::mr_presso(BetaOutcome = 'beta.outcome', BetaExposure = 'beta.exposure', SdOutcome = 'se.outcome',
                                          SdExposure = 'se.exposure', OUTLIERtest = TRUE, DISTORTIONtest = FALSE, data = x, NbDistribution = 1000,  SignifThreshold = 0.05)
        res_presso <- res_presso$`Main MR results`
        if(is.na(res_presso[nrow(res_presso),3])){
          res_presso <- res_presso[nrow(res_presso)-1,]
        } else {
          res_presso <- res_presso[nrow(res_presso),]
        }
        h_pval <- res_presso$`MR-PRESSO results`$`Global Test`$Pvalue
        res <- res %>% rbind(as_tibble(cbind(id.exposure=res$id.exposure[1],id.outcome=res$id.outcome[1],outcome=res$outcome[1],exposure=res$exposure[1],
                                             method='MRPRESSO',nsnp=res$nsnp[1],b=res_presso$`Causal Estimate`,se=res_presso$Sd,pval=res_presso$`P-value`)))
        
      }, silent = T
    )
  }
  if(is_empty(h_pval)){h_pval <- NA}
  res$F_stat <- MendelianRandomization::mr_ivw(x2[[1]])@Fstat
  res$Horizon_pval <- h_pval
  return(res)
})
mr_combined <- possibly(function(x){
  x2 <- dat_to_MRInput(x)
  res <- mr(x, method_list=ifelse(nrow(x)==1,'mr_wald_ratio','mr_ivw')) %>% as_tibble()
  res$heterogeneity_pval <- ifelse(nrow(x)==1, NA, mr_heterogeneity(x, method_list='mr_ivw')$Q_pval)
  pv <- mr_pleiotropy_test(x)$pval
  res$pleiotropy_pval <- ifelse(is_empty(pv), NA, pv)
  res$F_stat <- MendelianRandomization::mr_ivw(x2[[1]])@Fstat
  res$correct_causal_direction <- directionality_test(x)$correct_causal_direction
  return(res)
})
process_file <- possibly(function(e_data, o_file_name) {
  phenotype <- e_data$PROTEIN[1]
  e_data <- e_data %>% as_tibble %>% format_data(type = "exposure", chr_col = 'CHR', snp_col = 'SNP', beta_col = 'BETA', eaf_col = "MAF", se_col = 'SE', effect_allele_col = 'A2', other_allele_col = 'A1', pval_col = "P", phenotype_col = phenotype, samplesize_col = "N") %>% as_tibble
  o_data <- read_outcome_data(snps = e_data$SNP, filename = o_file_name, sep = "\t", chr_col = 'CHR', snp_col = "SNP", eaf_col = "MAF", beta_col = "BETA", se_col = "SE", pval_col = "P", effect_allele_col = "A1", other_allele_col = "A2", samplesize_col = "N") %>% as_tibble
  stopifnot(nrow(o_data)!=0)
  harmonise_data(exposure_dat = e_data, outcome_dat = o_data) %>% mr_combined %>% dplyr::mutate(exposure = phenotype, outcome = o_file_name) %>% select(-id.exposure, -id.outcome) %>% return()
})

full_mr_res <- NULL
factor_list <- file.path("/mnt/vol2/TMP_gwas_res/aging_traits/processed",paste0(aging_list,".tsv"))
for (i in 1:length(factor_list)) {
  trait <- basename(factor_list[i]) %>% gsub(".tsv","",.)
  print(paste0("Processing: ",trait))
  protein <- pQTL %>% filter(LOCATION == "cis") %>% split(.$PROTEIN)
  outcome_file <- factor_list[i]
  plan(multisession, workers = 20)
  full <- future_map(protein, ~ process_file(., outcome_file),.options = furrr_options(seed = TRUE)) %>% reduce(rbind) %>% arrange(pval) %>% mutate(P_BON = p.adjust(pval, method = "bonferroni"), P_FDR = p.adjust(pval, method = "BH"))
  full$outcome <- trait
  full_mr_res <- rbind(full_mr_res, full)
  plan(sequential)
}
full_mr_res %>% write_csv("result/step1_MR.csv")

full_mr_res <- NULL
#trait_assoc <- c('UMOD', 'HEXIM1', 'ACP1', 'EFEMP1', 'ATXN2L', 'FURIN', 'CELSR2', 'APOE', 'NT5C3A', 'IL1RN', 'KIT', 'SERPINF2', 'CEP170', 'SERPINA1', 'ARFIP1', 'NFATC1', 'TRDMT1', 'PARP1', 'TK1', 'TYMP', 'RPA2', 'COMMD1')
trait_assoc <- data.table::fread("result/step1_MR.csv") %>% filter(P_BON < 0.05) %>% pull(exposure) %>% unique
factor_list <- "/mnt/vol2/TMP_gwas_res/aging_traits/parent_lifespan.tsv"
for (i in 1:length(factor_list)) {
  trait <- basename(factor_list[i]) %>% gsub(".tsv","",.)
  print(paste0("Processing: ",trait))
  protein <- pQTL %>% filter(LOCATION == "cis", PROTEIN %in% trait_assoc) %>% split(.$PROTEIN)
  outcome_file <- factor_list[i]
  plan(multisession, workers = 20)
  full <- future_map(protein, ~ process_file(., outcome_file),.options = furrr_options(seed = TRUE)) %>% reduce(rbind) %>% arrange(pval) %>% mutate(P_BON = p.adjust(pval, method = "bonferroni"), P_FDR = p.adjust(pval, method = "BH"))
  full$outcome <- trait
  full_mr_res <- rbind(full_mr_res, full)
  plan(sequential)
}
full_mr_res %>% write_csv("result/step1_MR_parental.csv")
# Sens
sig_pro <- data.table::fread("result/step1_MR.csv") %>% filter(P_BON < 0.05)
pQTL <- pQTL %>% filter(PROTEIN %in% unique(sig_pro$exposure))

f <- function(x, y) {
  x <- pQTL %>% filter(PROTEIN == x)
  process_file(x, y)
}
trait <- sig_pro$outcome
trait_file <- file.path("/mnt/vol2/TMP_gwas_res/aging_traits/processed",paste0(trait,".tsv"))
protein <- sig_pro$exposure

plan(multisession, workers = 20)
full_mr_res <- future_pmap(list(protein, trait_file), ~ f(.x, .y), .options = furrr_options(seed = TRUE)) %>% reduce(rbind)
plan(sequential)

full_mr_res$outcome <- basename(full_mr_res$outcome) %>% gsub(".tsv","",.)
full_mr_res %>% write_csv("result/step1_MR_sens.csv")

# Validation
pQTL <- data.table::fread("/mnt/vol1/Resources/Finngen_pQTL/FinnGen_R10_pQTL_1e6_cis_trans_1Mb_0.001.tsv") %>% as_tibble
sig_pro <- data.table::fread("result/step1_MR.csv") %>% filter(P_BON < 0.05)
pQTL <- pQTL %>% filter(PROTEIN %in% unique(sig_pro$exposure))

f <- function(x, y) {
  x <- pQTL %>% filter(PROTEIN == x)
  process_file(x, y)
}
trait <- sig_pro$outcome
trait_file <- file.path("/mnt/vol2/TMP_gwas_res/aging_traits/processed",paste0(trait,".tsv"))
protein <- sig_pro$exposure

plan(multisession, workers = 20)
full_mr_res <- future_pmap(list(protein, trait_file), ~ f(.x, .y), .options = furrr_options(seed = TRUE)) %>% reduce(rbind)
plan(sequential)

full_mr_res$outcome <- gsub("/mnt/vol2/TMP_gwas_res/aging_traits/processed/","",full_mr_res$outcome)
full_mr_res$outcome <- gsub(".tsv","",full_mr_res$outcome)
full_mr_res$outcome[full_mr_res$outcome=="kdm_delta"] <- "KDM-BA acceleration"
full_mr_res$outcome[full_mr_res$outcome=="phenoage_delta"] <- "PhenoAge acceleration"
full_mr_res$outcome[full_mr_res$outcome=="ltl"] <- "Leukocyte telomere length"
full_mr_res$outcome[full_mr_res$outcome=="frailty_score"] <- "Frailty index"
full_mr_res$outcome[full_mr_res$outcome=="healthspan"] <- "Healthspan"
full_mr_res  %>% write_csv("result/step1_MR_finngen.csv")

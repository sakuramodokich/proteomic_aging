library(tidyverse)

# Extract data
map_ids <- c('p31','p21022','p21003_i0','p21000_i0','p22189','p20116_i0','p1558_i0','p21001_i0','p6138_i0','p54_i0',
             'p826_i0','p3426_i0','p30710_i0','p30000_i0','p30870_i0','p30750_i0','p30780_i0','p30650_i0','p30620_i0','p20003_i0','p1488_i0','p1498_i0','p30900_i0',
             'p1289_i0','p1299_i0','p1309_i0','p1319_i0','p1329_i0','p1339_i0','p1349_i0','p1369_i0','p1379_i0','p1389_i0','p1438_i0','p1448_i0','p1458_i0','p1468_i0')
cov <- data.table::fread("UKB_merged_dataset.tsv", select=c('eid', map_ids)) %>% as_tibble
cov$sex <- cov$p31
cov$age <- cov$p21022
cov$ethnic <- car::recode(cov$p21000_i0,'c(-1,-3)=NA;c(1,1001,1002,1003)=1;else=0')
cov$tdi <- cov$p22189
cov$qualification <- NA
cov$qualification[str_detect(cov$p6138_i0,"-7")] <- 0 # No
cov$qualification[str_detect(cov$p6138_i0,"2|3|4|5|6")] <- 1 # Other
cov$qualification[str_detect(cov$p6138_i0,"1")] <- 2 # High
cov$smoke_status <- car::recode(cov$p20116_i0, "-3=NA")
cov$alcohol_freq <- car::recode(cov$p1558_i0, "-3=NA;6=0;c(3,4,5)=1;c(1,2)=2")
cov$bmi <- cov$p21001_i0
cov$tg <- cov$p30870_i0
cov$hba1c <- cov$p30750_i0
cov$ass_center <- cov$p54_i0

food_item <- cov %>% select(starts_with(c("p1289_","p1299_","p1309_","p1319_","p1329_","p1339_","p1349_","p1369_","p1379_","p1389_","p1438_","p1448_","p1458_","p1468_"))) %>%
  mutate(across(everything(), ~case_when(.x %in% c(-1, -3) ~ NA, .x == -10 ~ 0.5, TRUE ~ as.numeric(.x))))
names(food_item) <- c('cooked_veg','raw_veg','fresh_fruit','dried_fruit','oily_fish','non_oily_fish','processed_meat','beef','lamb','pork',"bread","bread_type","cereal","cereal_type")
meat_var <- c('oily_fish', 'non_oily_fish', 'processed_meat', 'beef', 'lamb', 'pork')
food_item %>% mutate(
  across(all_of(meat_var), ~ case_when(.x == 1 ~ 0.5, .x == 2 ~ 1, .x == 4 ~ 5.5, .x == 5 ~ 7, TRUE ~ as.numeric(.x)))
) %>% mutate(
  all_veg = cooked_veg/2 + raw_veg,
  all_fruit = dried_fruit/2 + fresh_fruit,
  all_fish = oily_fish + non_oily_fish,
  unprocessed_meat = beef + lamb + pork,
  whole_grain = (bread * (bread_type %in% c(3, 4)) + cereal * (cereal_type %in% c(1, 4, 5)))/7,
  refined_grain = (bread * (bread_type %in% c(1, 2)) + cereal * (cereal_type %in% c(2, 3)))/7
) %>% mutate(
  all_veg = ifelse(all_veg>=3, 1, 0),
  all_fruit = ifelse(all_fruit>=3, 1, 0),
  all_fish = ifelse(all_fish>=2, 1, 0),
  processed_meat = ifelse(processed_meat<=1, 1, 0),
  unprocessed_meat = ifelse(unprocessed_meat<=1.5, 1, 0),
  whole_grain = ifelse(whole_grain>=3, 1, 0),
  refined_grain = ifelse(refined_grain<=1.5, 1, 0),
  diet_score = all_veg + all_fruit + processed_meat + unprocessed_meat + whole_grain + refined_grain + all_fish
) %>% .$diet_score -> cov$diet_score
rm(food_item, meat_var)

cov %>% select(eid, age, sex, tdi, ethnic, qualification, diet_score, smoke_status, alcohol_freq, bmi, ass_center) %>%
  mutate(
    qualification = factor(qualification),
    smoke_status = factor(smoke_status),
    alcohol_freq = factor(alcohol_freq)
  ) %>%
  write_csv("main_aging_cov.csv")

# GWAS
dat <- data.table::fread("main_aging_dataset.csv") %>% as_tibble
req <- c("p31","p22001","p22006","p22019","p22021","p22027","p21022","p22000",paste0("p22009_a",1:20))
cov <- data.table::fread("UKB_merged_dataset.tsv", select=c('eid', req)) %>% as_tibble
cov$age <- cov$p21022
cov$age_square <- cov$age ^ 2
cov$sex <- car::recode(cov$p31,"0=1;1=2")
cov$age_sex <- cov$sex * cov$age
cov$geno_batch <- cov$p22000
cov <- cov %>% filter(p31 == p22001, p22006 == 1, is.na(p22019), p22021 != 10, is.na(p22027))
cov <- cov %>% select(eid, all_of(paste0("p22009_a",1:20)), age, age_square, sex, age_sex, geno_batch)
names(cov)[2:21] <- paste0("pc",1:20)

dat <- dat %>% semi_join(cov, "eid")
cov <- cov %>% semi_join(dat, "eid")
dat <- dat %>% left_join(cov, "eid")
dat <- cbind(FID=dat$eid, IID=dat$eid, dat) %>% select(-eid)

pro <- data.table::fread("olink_protein_QC_knn_impute_INT.tsv") %>% pull(eid)
dat <- dat %>% filter(!IID %in% pro)
covs <- c('age', 'age_square', 'sex', 'age_sex', 'geno_batch', paste0("pc",1:20))

dat1 <- dat %>% drop_na(frailty_score) %>% select(FID, IID, frailty_score, all_of(covs))
dat1$frailty_score <-  RNOmni::RankNorm(dat1$frailty_score)

dat2 <- dat %>% drop_na(phenoage_delta) %>% select(FID, IID, phenoage_delta, all_of(covs))
dat2$phenoage_delta <-  RNOmni::RankNorm(dat2$phenoage_delta)

dat3 <- dat %>% drop_na(kdm_delta) %>% select(FID, IID, kdm_delta, all_of(covs))
dat3$kdm_delta <-  RNOmni::RankNorm(dat3$kdm_delta)

dat4 <- dat %>% drop_na(healthspan) %>% select(FID, IID, healthspan, all_of(covs))
dat4$healthspan <-  RNOmni::RankNorm(dat4$healthspan)

dat5 <- dat %>% drop_na(ltl) %>% select(FID, IID, ltl, all_of(covs))
dat5$ltl <-  RNOmni::RankNorm(dat5$ltl)

dat1 %>% write_delim("main_aging_gwas_frailty_score.txt", delim = " ")
dat2 %>% write_delim("main_aging_gwas_phenoage_delta.txt", delim = " ")
dat3 %>% write_delim("main_aging_gwas_kdm_delta.txt", delim = " ")
dat4 %>% write_delim("main_aging_gwas_healthspan_new.txt", delim = " ")
dat5 %>% write_delim("main_aging_gwas_ltl.txt", delim = " ")

library(tidyverse)
library(furrr)
library(BioAge)

req <- c(30620,30600,30610,30630,30640,30650,30710,30680,30690,30700,30720,30660,30730,30740,30750,
         30760,30770,30780,30790,30800,30810,30820,30830,30850,30840,30860,30870,30880,30670,30890,
         30180,30040,30070,30000,2306,22040,924,46,47,22191,31,4080,4079,20150,2247,6149,2453,53,34,52,
         2178,2080,1200,2050,1970,2020,1930,2188,2296,2463,2316,134,2335,6159,131354,21022,21003,3063)
dict <- data.table::fread("/mnt/vol2/UKB_dataset/asset/app79095_20231010100627.dataset.data_dictionary.csv") %>% filter(entity=="participant")
dict$field <- gsub("p","",dict$name) %>% str_split("[_]") %>% map_chr(1)
map_ids <- map(req, ~ dict$name[str_detect(dict$field, paste0("\\b",.,"\\b"))]) %>% unlist()
map_ids <- str_extract_all(map_ids, "^p\\d+$|^p\\d+_i0$|^p\\d+_i0_a\\d$") %>% unlist()
data_path <- "/mnt/vol2/UKB_dataset/asset/UKB_merged_dataset_23_10_recoded.tsv"
dat <- data.table::fread(data_path, select=c('eid', map_ids)) %>% as_tibble

dat$sbp <- rowMeans(dat[,c('p4080_i0_a0','p4080_i0_a1')], na.rm = T)
dat$sbp[is.nan(dat$sbp)] <- NA
dat$dbp <- rowMeans(dat[,c('p4079_i0_a0','p4079_i0_a1')], na.rm = T)
dat$dbp[is.nan(dat$dbp)] <- NA

dat$fev <- rowMeans(dat[,c('p3063_i0_a0','p3063_i0_a1','p3063_i0_a2')], na.rm = T)
dat$fev[is.nan(dat$fev)] <- NA

# KDM-BA
biomarkers_kdm = c("fev","sbp","totchol","hba1c","albumin","creat","lncrp","alp","bun")
biomarkers_phenoage = c("albumin_gL","lymph","mcv","glucose_mmol","rdw","creat_umol","lncrp","alp","wbc")

load("NHANES3.rda")
load("NHANES4.rda")
kdm = kdm_nhanes(biomarkers = biomarkers_kdm)

ukb <- dat %>% mutate(fev = fev * 1000, sbp = sbp, totchol = p30690_i0 * 38.67, hba1c = p30750_i0 * 0.0915 + 2.15, 
                      creat = p30700_i0 / 88.4, albumin = p30600_i0/10, lncrp = log(1 + p30710_i0/10), alp = p30610_i0, bun = p30670_i0 * 2.8, age = p21003_i0, sex = p31)
ukb_fem = kdm_calc(data = ukb %>% filter(sex==0), biomarkers_kdm, fit = kdm$fit$female, s_ba2 = kdm$fit$female$s_ba2)
ukb_male = kdm_calc(data = ukb %>% filter(sex==1), biomarkers_kdm, fit = kdm$fit$male, s_ba2 = kdm$fit$male$s_ba2)
ukb <- rbind(ukb_fem$data, ukb_male$data) %>% as_tibble %>% arrange(eid)

dat <- dat %>% left_join(ukb %>% select(eid, kdm, kdm_advance), "eid")
dat$kdm_delta <- dat$kdm_advance

# PhenoAge
phenoage = phenoage_nhanes(biomarkers_phenoage)

ukb <- dat %>% mutate(albumin_gL = p30600_i0, creat_umol = p30700_i0, glucose_mmol = p30740_i0, lncrp = log(1 + p30710_i0/10), lymph = p30180_i0,
                      mcv = p30040_i0, rdw = p30070_i0, alp = p30610_i0, wbc = p30000_i0, age = p21003_i0)
phenoage_ukb = phenoage_calc(data = ukb, biomarkers=biomarkers_phenoage, fit = phenoage$fit, orig = T) 
phenoage_ukb = phenoage_ukb$data
phenoage_ukb$phenoage[is.infinite(phenoage_ukb$phenoage)] <- NA
phenoage_ukb$phenoage_advance[is.infinite(phenoage_ukb$phenoage_advance)] <- NA

dat <- dat %>% left_join(phenoage_ukb %>% select(eid, phenoage, phenoage_advance), "eid")
dat$phenoage_delta <- dat$phenoage_advance

# Frailty
dat$birth <- as.Date(paste(dat$p34, dat$p52, "15", sep = "-"), "%Y-%m-%d")
if(F){
  self_non_ca <- data.table::fread(data_path, select=c('eid', 'p20002_i0', paste0('p20008_i0_a',0:33))) %>% as_tibble
  self_non_ca$p20002_i0[self_non_ca$p20002_i0 == ""] <- NA
  lst <- str_split(self_non_ca$p20002_i0, "[|]")
  max_length <- 34
  pad <- function(x, max_length, pad = NA){
    length(x) <- max_length
    x
  }
  padded_list <- map(lst, ~ pad(.x, max_length, pad = NA))
  tibble_df <- future_map_dfr(padded_list, as_tibble, .options = furrr_options(seed = T))
  dd <- self_non_ca %>% select(-p20002_i0) %>% pivot_longer(cols = starts_with("p20008_"), names_to = "key", values_to = "value")
  dd <- cbind(dd, tibble_df)
  dd <- dd %>% drop_na()
  dd <- dd[,-2]
  names(dd) <- c("eid","epistart","code")
  dd %>% write_csv("Revise/non_cancer_self_reported.csv")
}
self_non_ca <- data.table::fread("Revise/non_cancer_self_reported.csv") %>% as_tibble
self_non_ca$epistart <- car::recode(self_non_ca$epistart,"c(-1,-3)=NA")
self_non_ca <- self_non_ca %>% drop_na(epistart)

self_non_ca$year_birth <- dat$p34[match(self_non_ca$eid,dat$eid)]
self_non_ca$mo_birth <- dat$p52[match(self_non_ca$eid,dat$eid)]
self_non_ca$year_birth <- self_non_ca$year_birth + round((0.1 + (self_non_ca$mo_birth - 1) * 0.8/11),1)
self_non_ca$epistart <- self_non_ca$epistart - self_non_ca$year_birth

map_self <- function(data, map_code){
  map_set <- self_non_ca[str_detect(self_non_ca$code, map_code),] %>% arrange(eid, epistart) %>% filter(!duplicated(eid))
  data$epistart <- NA
  data$epistart[match(map_set$eid, data$eid)] <- map_set$epistart
  data$epistart <- as.numeric((data$epistart - data$p21003_i0)<=0)
  data$epistart[is.na(data$epistart)] <- 0
  return(data$epistart)
}
dat$glaucoma <- map_self(dat, "1277")
dat$cataracts <- map_self(dat, "1278")
dat$migraine <- map_self(dat, "1265")
dat$anxiety_panic <- map_self(dat, "1287")
dat$diabetes <- map_self(dat, "1220")
dat$myocardial_infarction <- map_self(dat, "1075")
dat$angina <- map_self(dat, "1074")
dat$stroke <- map_self(dat, "1081")
dat$hypertension <- map_self(dat, "1065")
dat$hypothyroidism <- map_self(dat, "1226")
dat$deep_vein_thrombosis <- map_self(dat, "1094")
dat$high_cholesterol <- map_self(dat, "1473")
dat$pneumonia <- map_self(dat, "1398")
dat$chronic_bronchitis <- map_self(dat, "1113")
dat$asthma <- map_self(dat, "1111")
dat$rheumatoid_arthritis <- map_self(dat, "1464")
dat$osteoarthritis <- map_self(dat, "1465")
dat$gout <- map_self(dat, "1466")
dat$osteoporosis <- map_self(dat, "1309")
dat$hayfever_allergic_rhinitis_eczema <- map_self(dat, "1387|1452")
dat$psoriasis <- map_self(dat, "1453")
dat$sciatica <- map_self(dat, "1476")
dat$gastric_reflux <- map_self(dat, "1138")
dat$hiatus_hernia <- map_self(dat, "1474")
dat$gall_stones <- map_self(dat, "1162")
dat$diverticulitis <- map_self(dat, "1458")

dat$cancer_doctor <- car::recode(dat$p2453_i0, "c(-1,-3)=NA")
dat$hearing_difficulty <- car::recode(dat$p2247_i0, "c(-1,-3)=NA;99=1")
dat$dental_problems <- car::recode(dat$p6149_i0, "-3=NA;c('',-7)=0;NA=NA;else=1")
dat$self_rated_health <- car::recode(dat$p2178_i0, "c(-1,-3)=NA;1=0;2=0.25;3=0.5;4=1")
dat$fatigue <- car::recode(dat$p2080_i0, "c(-1,-3)=NA;1=0;2=0.25;3=0.5;4=1")
dat$sleep_insomnia <- car::recode(dat$p1200_i0, "-3=NA;1=0;2=0.5;3=1")
dat$depressed_feeling <- car::recode(dat$p2050_i0, "c(-1,-3)=NA;1=0;2=0.5;3=0.75;4=1")
dat$nervous_feeling <- car::recode(dat$p1970_i0, "c(-1,-3)=NA")
dat$loneliness_feeling <- car::recode(dat$p2020_i0, "c(-1,-3)=NA")
dat$misery_feeling <- car::recode(dat$p1930_i0, "c(-1,-3)=NA")
dat$longstanding_ill <- car::recode(dat$p2188_i0, "c(-1,-3)=NA")
dat$fall_lastyear <- car::recode(dat$p2296_i0, "-3=NA;1=0;2=0.5;3=1")
dat$fractures_lastfy <- car::recode(dat$p2463_i0, "c(-1,-3)=NA")
dat$wheeze_lastyear <- car::recode(dat$p2316_i0, "c(-1,-3)=NA")
dat$multiple_cancer <- car::recode(dat$p134_i0, "0=0;NA=NA;else=1")
dat$chest_pain <- car::recode(dat$p2335_i0, "c(-1,-3)=NA")

dat$p6159_i0 <- car::recode(dat$p6159_i0, "-3=NA")
dat$head_neck_pain <- as.numeric(dat$p6159_i0 == 1 | dat$p6159_i0 == 3)
dat$back_pain <- as.numeric(dat$p6159_i0 == 4)
dat$stomach_pain <- as.numeric(dat$p6159_i0 == 5)
dat$hip_pain <- as.numeric(dat$p6159_i0 == 6)
dat$knee_pain <- as.numeric(dat$p6159_i0 == 7)
dat$whole_body_pain <- as.numeric(dat$p6159_i0 == 8)
dat$facial_pain <- as.numeric(dat$p6159_i0 == 2)

fi_item <- c('glaucoma', 'cataracts', 'migraine', 'anxiety_panic', 'diabetes',
             'myocardial_infarction', 'angina', 'stroke', 'hypertension', 'hypothyroidism',
             'deep_vein_thrombosis', 'high_cholesterol', 'pneumonia', 'chronic_bronchitis',
             'asthma', 'rheumatoid_arthritis', 'osteoarthritis', 'gout', 'osteoporosis',
             'hayfever_allergic_rhinitis_eczema', 'psoriasis', 'sciatica', 'gastric_reflux',
             'hiatus_hernia', 'gall_stones', 'diverticulitis', 'cancer_doctor', 'hearing_difficulty', 
             'dental_problems', 'self_rated_health', 'fatigue', 'sleep_insomnia', 'depressed_feeling', 'nervous_feeling',
             'loneliness_feeling', 'misery_feeling', 'longstanding_ill', 'fall_lastyear', 'fractures_lastfy',
             'wheeze_lastyear', 'multiple_cancer', 'chest_pain', 'head_neck_pain', 'back_pain', 'stomach_pain',
             'hip_pain', 'knee_pain', 'whole_body_pain', 'facial_pain')
dat$frailty_item_missing <- matrixStats::rowCounts(dat %>% select(all_of(fi_item)) %>% as.matrix, value = NA_integer_)
imp <- mice::mice(dat %>% filter(frailty_item_missing >= 0, frailty_item_missing <= 9) %>% select(eid, all_of(fi_item)) %>% mutate(across(-eid, factor)),
                  maxit = 5, m = 5, seed=123)
result <- complete(imp)
dat <- dat %>% select(-all_of(fi_item)) %>% left_join(result, "eid")
dat <- dat %>% 
  mutate_if(is.factor, as.numeric) %>%
  mutate(frailty_score = (glaucoma + cataracts + migraine + anxiety_panic + diabetes +
                            myocardial_infarction + angina + stroke + hypertension + hypothyroidism +
                            deep_vein_thrombosis + high_cholesterol + pneumonia + chronic_bronchitis +
                            asthma + rheumatoid_arthritis + osteoarthritis + gout + osteoporosis +
                            hayfever_allergic_rhinitis_eczema + psoriasis + sciatica + gastric_reflux +
                            hiatus_hernia + gall_stones + diverticulitis + cancer_doctor + hearing_difficulty + 
                            dental_problems + self_rated_health + fatigue + sleep_insomnia + depressed_feeling + nervous_feeling +
                            loneliness_feeling + misery_feeling + longstanding_ill + fall_lastyear + fractures_lastfy +
                            wheeze_lastyear + multiple_cancer + chest_pain + head_neck_pain + back_pain + stomach_pain +
                            hip_pain + knee_pain + whole_body_pain + facial_pain)/49)

# Telomere length
dat$ltl <- dat$p22191_i0

# Healthspan
icd10 <- data.table::fread("~/vol2/UKB_dataset/health_record/inpatient_diagnosis.csv")

death <- data.table::fread("~/vol2/UKB_dataset/health_record/death.csv") %>% arrange(eid)
cancer <- data.table::fread("Revise/cancer.csv")
cancer <- data.table::fread(data_path, select=names(cancer)) %>% as_tibble
cancer <- car::recode(cancer,"c(-1,-3)=NA")

cancer$cancer <- matrixStats::rowMins(as.matrix(cancer[,-1]),na.rm = T)
cancer$cancer[is.infinite(cancer$cancer)] <- NA
cancer$cancer <- round(cancer$cancer,0)

map_outcome <- function(data, map_code){
  map_set <- icd10[str_detect(icd10$code, map_code),] %>% arrange(eid, epistart) %>% filter(!duplicated(eid))
  data$epistart <- data.table::as.IDate(NA)
  data$epistart[match(map_set$eid, data$eid)] <- map_set$epistart
  data$epistart <- round(as.numeric(difftime(data$epistart, data$birth, units = "days"))/365.25,1)
  #data$epistart[which(data$epistart - data$p21003_i0 > 0)] <- NA
  return(data$epistart)
}
dat$CHF <- map_outcome(dat, "I50")
dat$COPD <- map_outcome(dat, "J44")
dat$MI <- map_outcome(dat, "I2[1-5]")
dat$DEMENTIA <- map_outcome(dat, "F0[1-5]")
dat$DIABETES <- map_outcome(dat, "E1[0-4]")
dat$STROKE <- map_outcome(dat, "I6[0-4]")
dat$CANCER <- map_outcome(dat, "C[0-9][0-9](?!44)")

map_self <- function(data, map_code){
  map_set <- self_non_ca[str_detect(self_non_ca$code, map_code),] %>% arrange(eid, epistart) %>% filter(!duplicated(eid))
  data$epistart <- NA
  data$epistart[match(map_set$eid, data$eid)] <- map_set$epistart
  #data$epistart[data$epistart - data$p21003_i0 > 0] <- NA
  return(data$epistart)
}
dat$CHF_sr <- map_self(dat, "1076")
dat$COPD_sr <- map_self(dat, "1112")
dat$MI_sr <- map_self(dat, "1075")
dat$DEMENTIA_sr <- map_self(dat, "1263|1258|1259|1260|1261|1262")
dat$DIABETES_sr <- map_self(dat, "1220|1221|1222|1223|1521")
dat$STROKE_sr <- map_self(dat, "1081|1086|1491|1583")

dat$death <- data.table::as.IDate(NA)
dat$death[match(death$eid,dat$eid)] <- death$date_of_death
dat$death <- round(as.numeric(difftime(dat$death, dat$birth, units = "days"))/365.25,0)

dat$cancer_regi <- cancer$cancer

dat$healthspan <- NA
disease <- dat %>% select(CHF,COPD,MI,DEMENTIA,DIABETES,STROKE,CANCER,CHF_sr,COPD_sr,MI_sr,DEMENTIA_sr,DIABETES_sr,STROKE_sr, death,cancer_regi)
disease <- disease %>%
  mutate(
    CHF = pmin(CHF, CHF_sr, na.rm = T),
    COPD = pmin(COPD, COPD_sr, na.rm = T),
    MI = pmin(MI, MI_sr, na.rm = T),
    DEMENTIA = pmin(DEMENTIA, DEMENTIA_sr, na.rm = T),
    DIABETES = pmin(DIABETES, DIABETES_sr, na.rm = T),
    STROKE = pmin(STROKE, STROKE_sr, na.rm = T),
    CANCER = pmin(CANCER, cancer_regi, na.rm = T)
  )
disease$sum <- matrixStats::rowMins(as.matrix(disease), na.rm = T)
disease$sum[is.infinite(disease$sum)] <- NA
dat$healthspan <- disease$sum
dat$year_birth <- dat$p21003_i0 + round((0.1 + (dat$p52 - 1) * 0.8/11),1)

dat$healthspan[is.na(dat$healthspan)] <- dat$year_birth[is.na(dat$healthspan)]

dat <- dat %>% select(eid, healthspan)
dat %>% write_csv("main_aging_dataset_new.csv")

### Forward Mendelian randomization ###
expdata <- data.table::fread(file = "input_raw_data.txt")
expdata_fil <-subset(expdata,P<1e-05)
expdata_fil$beta <- log(expdata_fil$OR)
write.csv(expdata_fil, file="fill_data.csv")
expdata_unclumped <- read_exposure_data(filename = "fill_data.csv",
                                        sep = ",",
                                        snp_col = "SNP",
                                        beta_col = "beta",
                                        se_col = "se",
                                        eaf_col = "eaf",
                                        effect_allele_col = "effect_allele",
                                        other_allele_col = "other_allele",
                                        clump = F) 
expdata_clumped <-clump_data(expdata_unclumped,
                             clump_kb = 10000,
                             clump_r2 = 0.0001,)
save(expdata_clumped,file = "clumped_data.Rdata")
load("tempdata/MDD_expdata_clumped.Rdata")
outcomedata <- extract_outcome_data(snps=expdata_clumped$SNP,
                            outcomes='outcome_id', 
                            proxies = FALSE,
                            maf_threshold = 0.01,
                            access_token = NULL)
dat <-harmonise_data(exposure_dat = expdata_clumped,
                     outcome_dat = outcomedata)
dat_MRres <- mr(dat)
dat_MRres_addOR <- generate_odds_ratios(mr_res = dat_MRres)
mr_scatter_plot(mr_results = dat_MRres,dat)
heterogeneity <- mr_heterogeneity(dat)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
pleiotropy <- mr_pleiotropy_test(dat)
leaveoneout_results <- mr_leaveoneout(dat)
mr_leaveoneout_plot(leaveoneout_results)
mr_presso <- run_mr_presso(dat,NbDistribution = 3000)
mr_presso

### Reverse Mendelian randomization ###
library(TwoSampleMR)
expdata <- extract_instrumenoutcome(
  outcomes='Exposure_id',
  p1=1e-5,
  clump=TRUE,
  r2=0.0001,
  kb=10000,
  access_token = NULL)
save(expdata,file = "expdata_clumped.Rdata")
load("expdata_clumped.Rdata")
outcome <-data.table::fread("outcome_raw_data.txt",header = T)
head(outcome)
outcome <-merge(expdata,outcome,by.x = "SNP",by.y = "SNP")
write.csv(outcome, file = "outcome.csv")
outcome_clean <- read_outcome_data(filename = "outcome.csv",
                                   sep = ",",
                                   snps = expdata$SNP,
                                   snp_col = "SNP",
                                   beta_col = "beta",
                                   se_col = "se",
                                   effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele",
                                   pval_col = "Pvalue")
dat <-harmonise_data(exposure_dat = expdata,
                     outcome_dat = outcome_clean)
dat_MRres <- mr(dat)
dat_MRres_addOR <- generate_odds_ratios(mr_res = dat_MRres)
mr_scatter_plot(mr_resuloutcome = dat_MRres,dat)
heterogeneity <- mr_heterogeneity(dat)
mr_funnel_plot(singlesnp_resuloutcome = mr_singlesnp(dat))
pleiotropy <- mr_pleiotropy_test(dat)
leaveoneout_results <- mr_leaveoneout(dat)
mr_leaveoneout_plot(leaveoneout_results)
mr_presso <- run_mr_presso(dat,NbDistribution = 3000)
mr_presso
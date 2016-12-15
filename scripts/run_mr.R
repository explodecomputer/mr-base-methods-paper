library(TwoSampleMR)
library(MRInstruments)
library(dplyr)

# 1. Get list of outcomes to include
# 2. Define instruments
# 3. Extract outcome data
# 4. Perform MR
# 5. Make plots



# Get all available outcomes
# toggle_dev("test")
ao <- available_outcomes()


# Select most reliable studies for diseases and risk factors
# Exclude lipids
outcomes <- filter(ao, 
		category %in% c("Disease", "Risk factor"),
		(grepl("Mixed", population) | grepl("European", population)),
		mr == 1,
		sample_size > 1000,
		nsnp > 95000,
		access %in% c("immunobase_users", "public")
	) %>%
	arrange(desc(year), desc(nsnp)) %>%
	distinct(trait, .keep_all=TRUE) %>%
	filter(! id %in% c(299:302))



# Get Lp(a) instrument
lpa <- read_exposure_data("../data/lpa_inst.txt")



# Get trigs (302) and LDL (300) instruments, Use GLGC list of robust associations, Extract those from the GLGC summary data, Keep p < 5e-8 and clump
# glgc_snps <- scan("../data/glgc_snps.txt", what="character") %>%
# 	extract_outcome_data(snps=., outcomes = c(300, 302)) %>%
# 	convert_outcome_to_exposure() %>%
# 	filter(pval.exposure < 5e-8) %>%
# 	clump_data()

data(gwas_catalog)

ldl_inst <- gwas_catalog %>%
	filter(Phenotype == "LDL cholesterol (unit increase)", Author == "Willer CJ", Year == 2013) %>%
	format_data()

trig_inst <- gwas_catalog %>%
	filter(Phenotype == "Triglycerides (mg/dL increase)", Author == "Willer CJ", Year == 2013) %>%
	format_data()





# Exposure data
exposure_dat <- combine_data(list(ldl_inst, trig_inst, lpa))
exposure_dat$exposure[grepl("LDL", exposure_dat$exposure)] <- "LDL cholesterol"
exposure_dat$exposure[grepl("Trig", exposure_dat$exposure)] <- "Triglycerides"
exposure_dat$other_allele.exposure[exposure_dat$SNP == "rs2072183"] <- "G"



# Extract outcome data
outcome_dat <- extract_outcome_data(exposure_dat$SNP, outcomes$id)


# Set units to SD
outcome_dat <- subset(outcome_dat, id.outcome %in% outcomes$id)
outcomes <- subset(outcomes, id %in% outcome_dat$id.outcome)
index <- ! grepl("SD", outcomes$unit, ignore.case=TRUE) & ! grepl("log odds", outcomes$unit, ignore.case=TRUE) & ! is.na(outcomes$sd)
temp <- outcomes[index,] %>% select(id, sd)
outcome_dat <- merge(outcome_dat, temp, by.x="id.outcome", by.y="id", all.x=TRUE)
index <- !is.na(outcome_dat$sd)
outcome_dat$beta.outcome[index] <- outcome_dat$beta.outcome[index] / outcome_dat$sd[index]
outcome_dat$se.outcome[index] <- outcome_dat$se.outcome[index] / outcome_dat$sd[index]


# Flip Asthma - the database needs to be updated
outcome_dat$beta.outcome[grepl("Asthma", outcome_dat$outcome)] <- outcome_dat$beta.outcome[grepl("Asthma", outcome_dat$outcome)] * -1

# Harmonise data
dat1 <- harmonise_data(exposure_dat, subset(outcome_dat, ! id.outcome %in% c(23,784)), action=2)
dat2 <- harmonise_data(exposure_dat, subset(outcome_dat, id.outcome %in% c(23,784)), action=1)
dat <- rbind(dat1, dat2)


# Perform MR
res <- mr(dat)


# Flip results in MR to reflect influence of intervention
res$b <- res$b * -1


# Save results
save(outcomes, dat, res, exposure_dat, outcome_dat, file="../results/lpa_ldl_trigs.RData")


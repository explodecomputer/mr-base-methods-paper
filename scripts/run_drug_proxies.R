library(TwoSampleMR)
library(dplyr)
toggle_dev("test")


# Get LDL drug proxies
# rs10790162
# rs11591147
ldl_drug_proxies <- extract_outcome_data(c("rs12916", "rs11591147", "rs2073547"), 300)
ldl_drug_proxies$gene <- c("HMGCR", "PCSK9", "NPC1L1")
ldl_drug_proxies$outcome <- c("HMGCR inhibitor", "PCSK9 inhibitor", "NPC1L1 inhibitor")


# Get trig drug proxies
trig_drug_proxies <- extract_outcome_data("rs10790162", 302)
trig_drug_proxies$gene <- "APOC3"
trig_drug_proxies$outcome <- "APOC3 inhibitor"


# Get LPA drug proxies
# P value from Ye PMID: 24089516, standard error not published, calculated using z val
lpa_drug_proxies <- read_exposure_data("../data/lpa_inst.txt")
lpa_drug_proxies$gene.exposure <- "LPA"
lpa_drug_proxies$exposure <- "Lp(a) inhibitor"


# Create exposure data
exposure_dat <- bind_rows(
	convert_outcome_to_exposure(ldl_drug_proxies), 
	convert_outcome_to_exposure(trig_drug_proxies), 
	lpa_drug_proxies
)
exposure_dat$units.exposure <- "SD"


# Select most reliable studies for diseases and risk factors
# Exclude lipids
ao <- available_outcomes()
outcomes1 <- filter(ao, 
		category %in% c("Disease", "Risk factor"),
		(grepl("Mixed", population) | grepl("European", population)),
		mr == 1,
		sample_size > 1000,
		nsnp > 95000
	) %>%
	arrange(desc(year), desc(nsnp)) %>%
	distinct(trait, .keep_all=TRUE)
outcomes2 <- filter(ao, id == 24)
outcomes <- rbind(outcomes1, outcomes2)


# Extract outcome data
outcome_dat <- extract_outcome_data(exposure_dat$SNP, outcomes$id)


# Harmonise data
dat1 <- harmonise_data(exposure_dat, subset(outcome_dat, id.outcome != 784), action=2)
dat2 <- harmonise_data(exposure_dat, subset(outcome_dat, id.outcome == 784), action=1)
dat <- rbind(dat1, dat2)


# Perform MR
res <- mr(dat)


# Flip results in MR to reflect influence of intervention
res$b <- res$b * -1


# Save results
save(outcomes, dat, res, exposure_dat, outcome_dat, file="../results/drug_proxies.RData")



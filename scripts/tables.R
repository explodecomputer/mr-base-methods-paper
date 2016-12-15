library(TwoSampleMR)
library(dplyr)
library(RMySQL)
library(googlesheets)
# toggle_dev("test")

## Studies

ao <- available_outcomes()
sheet1 <- filter(ao, pmid != 0, access %in% c("public", "Immunobase")) %>%
	group_by(pmid) %>%
	summarise(
		"Pubmed ID" = first(pmid),
		Phenotypes = paste(unique(trait), collapse=", "),
		"First author" = first(author),
		Year = first(year),
		Consortium = ifelse(is.na(first(consortium)), "", first(consortium)),
		"No. GWAS analyses" = n(),
		"No. cases" = first(ncase),
		"No. controls" = first(ncontrol),
		"Sample size" = first(sample_size),
		Category = first(category),
		Subcategory = first(subcategory),
		Access = first(access)
	) %>%
	arrange(Access, Phenotypes) %>%
	select(-pmid)



# Drug prediction

load("../results/drug_proxies.RData")
load("../results/drug_prediction.RData")
# res <- select(res, test, exposure, outcome, b, se)
lab <- do.call(rbind, strsplit(as.character(res$labels), split="\n"))
res$pval <- pnorm(res$b/res$se)
res$target <- gsub(" variant", "", lab[,2])
res$instrument <- lab[,1]
res$pathway <- gsub("Target: ", "", lab[,3])
res$ref <- NA

res_mr <- subset(res, test=="Mendelian randomization")
res_rct <- subset(res, test!="Mendelian randomization")

inst <- subset(exposure_dat, select=c(SNP, exposure))
res_mr <- merge(res_mr, inst, by="exposure")
res_mr$instrument <- res_mr$SNP
res_mr <- subset(res_mr, select=-c(SNP))

res_mr$ref[res_mr$exposure == "APOC3 inhibitor" & res_mr$outcome == "Type 2 diabetes"] <- "GLGC 2013, DIAGRAM 2014"
res_mr$ref[res_mr$exposure == "HMGCR inhibitor" & res_mr$outcome == "Type 2 diabetes"] <- "GLGC 2013, DIAGRAM 2014"
res_mr$ref[res_mr$exposure == "Lp(a) inhibitor" & res_mr$outcome == "Type 2 diabetes"] <- "Qi et al 2012, DIAGRAM 2014"
res_mr$ref[res_mr$exposure == "NPC1L1 inhibitor" & res_mr$outcome == "Type 2 diabetes"] <- "GLGC 2013, DIAGRAM 2014"
res_mr$ref[res_mr$exposure == "PCSK9 inhibitor" & res_mr$outcome == "Type 2 diabetes"] <- "GLGC 2013, DIAGRAM 2012"
res_mr$ref[res_mr$exposure == "APOC3 inhibitor" & res_mr$outcome == "Coronary heart disease"] <- "GLGC 2013, Nikpay et al 2015"
res_mr$ref[res_mr$exposure == "HMGCR inhibitor" & res_mr$outcome == "Coronary heart disease"] <- "GLGC 2013, Nikpay et al 2015"
res_mr$ref[res_mr$exposure == "Lp(a) inhibitor" & res_mr$outcome == "Coronary heart disease"] <- "Qi et al 2012, Nikpay et al 2015"
res_mr$ref[res_mr$exposure == "NPC1L1 inhibitor" & res_mr$outcome == "Coronary heart disease"] <- "GLGC 2013, Nikpay et al 2015"
res_mr$ref[res_mr$exposure == "PCSK9 inhibitor" & res_mr$outcome == "Coronary heart disease"] <- "GLGC 2013, Nikpay et al 2015"


res_rct$ref[res_rct$exposure == "HMGCR inhibitor" & res_rct$outcome == "Type 2 diabetes"] <- "Swerdlow et al 2015"
res_rct$ref[res_rct$exposure == "NPC1L1 inhibitor" & res_rct$outcome == "Type 2 diabetes"] <- "Cannon et al 2015"
res_rct$ref[res_rct$exposure == "PCSK9 inhibitor" & res_rct$outcome == "Type 2 diabetes"] <- "Sabatine et al 2015"
res_rct$ref[res_rct$exposure == "HMGCR inhibitor" & res_rct$outcome == "Coronary heart disease"] <- "Ference et al 2012"
res_rct$ref[res_rct$exposure == "NPC1L1 inhibitor" & res_rct$outcome == "Coronary heart disease"] <- "Cannon et al 2015"
res_rct$ref[res_rct$exposure == "PCSK9 inhibitor" & res_rct$outcome == "Coronary heart disease"] <- "Sabatine et al 2015"

res_rct <- rbind(res_rct, res_mr[res_mr$exposure == "APOC3 inhibitor" & res_mr$outcome == "Type 2 diabetes", ])
res_rct$test[nrow(res_rct)] <- "Randomized controlled trial"
res_rct$pval[nrow(res_rct)] <- NA
res_rct$b[nrow(res_rct)] <- NA
res_rct$se[nrow(res_rct)] <- NA
res_rct$ref[nrow(res_rct)] <- NA

res_rct <- rbind(res_rct, res_mr[res_mr$exposure == "Lp(a) inhibitor" & res_mr$outcome == "Type 2 diabetes", ])
res_rct$test[nrow(res_rct)] <- "Randomized controlled trial"
res_rct$pval[nrow(res_rct)] <- NA
res_rct$b[nrow(res_rct)] <- NA
res_rct$se[nrow(res_rct)] <- NA
res_rct$ref[nrow(res_rct)] <- NA

res_rct <- rbind(res_rct, res_mr[res_mr$exposure == "APOC3 inhibitor" & res_mr$outcome == "Coronary heart disease", ])
res_rct$test[nrow(res_rct)] <- "Randomized controlled trial"
res_rct$pval[nrow(res_rct)] <- NA
res_rct$b[nrow(res_rct)] <- NA
res_rct$se[nrow(res_rct)] <- NA
res_rct$ref[nrow(res_rct)] <- NA

res_rct <- rbind(res_rct, res_mr[res_mr$exposure == "Lp(a) inhibitor" & res_mr$outcome == "Coronary heart disease", ])
res_rct$test[nrow(res_rct)] <- "Randomized controlled trial"
res_rct$pval[nrow(res_rct)] <- NA
res_rct$b[nrow(res_rct)] <- NA
res_rct$se[nrow(res_rct)] <- NA
res_rct$ref[nrow(res_rct)] <- NA


sheet3 <- rbind(res_rct, res_mr)
sheet3 <- subset(sheet3, select=c(exposure, outcome, target, instrument, pathway, test, b, se, pval, ref))
sheet3 <- sheet3[order(sheet3$outcome, sheet3$target),]



# Instruments

ucsc_get_position <- function(snp)
{
	snp <- paste(snp, collapse="', '")
	require(RMySQL)
	message("Connecting to UCSC MySQL database")
	mydb <- dbConnect(MySQL(), user="genome", dbname="hg19", host="genome-mysql.cse.ucsc.edu")

	query <- paste0(
		"SELECT * from snp144 where name in ('", snp, "');"
	)
	message(query)
	out <- dbSendQuery(mydb, query)
	d <- fetch(out, n=-1)
	# dbClearResult(dbListResults(mydb)[[1]])
	dbDisconnect(mydb)
	return(d)
}


load("../results/lpa_ldl_trigs.RData")

pos <- ucsc_get_position(unique(exposure_dat$SNP)) %>%
	select(name, chrom, chromStart)

exposure_dat <- merge(exposure_dat, pos, by.x="SNP", by.y="name")
names(exposure_dat) <- gsub(".exposure", "", names(exposure_dat))

sheet4 <- select(exposure_dat,
		exposure,
		SNP,
		gene,
		chrom,
		chromStart,
		effect_allele,
		other_allele,
		eaf,
		beta,
		se,
		pval
	) %>%
	mutate(chrom = sapply(strsplit(chrom, split="_"), function(x) as.numeric(gsub("chr", "", x[[1]]))
		)) %>%
	arrange(exposure, chrom, chromStart)



# Main results

load("../results/lpa_ldl_trigs.RData")
outs <- select(outcomes, id, trait, pmid, year, author, consortium, subcategory, ncase, ncontrol, sample_size, unit, sd, sex, population)

res <- filter(res, (nsnp != 1 & method %in% c("MR Egger", "Inverse variance weighted", "Weighted median")) | (nsnp == 1 & method == "Wald ratio"))

het <- mr_heterogeneity(dat, method_list=c("mr_egger_regression", "mr_ivw")) %>% 
	select(id.exposure, id.outcome, method, Q, Q_df, Q_pval) 
plei <- mr_pleiotropy_test(dat) %>% transmute(id.exposure=id.exposure, id.outcome=id.outcome, egger_intercept=egger_intercept, egger_intercept_se=se, egger_intercept_pval=pval)
plei$method <- "MR Egger"

tab <- full_join(het, plei)
tab <- full_join(res, tab)

sheet5 <- merge(tab, outs, by.x="id.outcome", by.y="id") %>%
	mutate(outcome=trait) %>%
	select(-id.outcome, -id.exposure, -trait) %>%
	arrange(exposure, subcategory, outcome, method)




# Drug proxy results

load("../results/drug_proxies.RData")
outs <- select(outcomes, id, trait, pmid, year, author, consortium, subcategory, ncase, ncontrol, sample_size, unit, sex, population)

sheet6 <- merge(res, outs, by.x="id.outcome", by.y="id") %>%
	mutate(outcome=trait) %>%
	select(-id.outcome, -id.exposure, -trait, -method) %>%
	arrange(exposure, subcategory, outcome)




# Write out

write.csv(sheet1, file="../tables/study.csv")
write.csv(sheet3, file="../tables/drug_prediction.csv")
write.csv(sheet4, file="../tables/instruments.csv")
write.csv(sheet5, file="../tables/mr_results.csv")
write.csv(sheet6, file="../tables/drug_proxy_results.csv")



# Upload

# gs <- gs_title("MR-Base methods paper supplementary tables")
# gs_delete(gs)
gs <- gs_new(title="MR-Base methods paper supplementary tables")
gs_ws_new(gs, ws_title = "Supplementary table 1", input = sheet1)
gs_ws_new(gs, ws_title = "Supplementary table 3", input = sheet3)
gs_ws_new(gs, ws_title = "Supplementary table 4", input = sheet4)
gs_ws_new(gs, ws_title = "Supplementary table 5", input = sheet5)
# gs_ws_new(gs, ws_title = "Supplementary table 6", input = sheet6)


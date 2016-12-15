library(dplyr)
load("../results/lpa_ldl_trigs.RData")

a <- subset(res, exposure == "LDL cholesterol" & grepl("Type 2", outcome))
a$lor <- exp(a$b)
a$lo <- exp(a$b - 1.96 * a$se)
a$up <- exp(a$b + 1.96 * a$se)


table(exposure_dat$exposure)


table(subset(res, method == "MR Egger")$pval < 0.05)
table(subset(res, method == "Weighted median")$pval < 0.05)
table(subset(res, method == "Inverse variance weighted")$pval < 0.05)

sum(!is.na(subset(res, method == "MR Egger")$pval))
sum(!is.na(subset(res, method == "Weighted median")$pval))
sum(!is.na(subset(res, method == "Inverse variance weighted")$pval))

sum(subset(res, method == "Inverse variance weighted")$pval < 0.05) %>% binom.test(n=294, p=0.05)

sum(subset(res, method %in% c("Wald ratio", "Inverse variance weighted"))$pval < 0.05) %>% binom.test(n=length(subset(res, method %in% c("Wald ratio", "Inverse variance weighted"))$pval < 0.05), p=0.05)


sum(subset(res, method == "Weighted median")$pval < 0.05) %>% binom.test(n=294, p=0.05)

sum(subset(res, method == "MR Egger")$pval < 0.05) %>% binom.test(n=294, p=0.05)



table(subset(res, exposure == "LDL cholesterol" & method == "Inverse variance weighted")$pval < 0.05) %>% prop.test(p=c(0.95))

table(subset(res, exposure == "Triglycerides" & method == "Inverse variance weighted")$pval < 0.05) %>% prop.test(p=c(0.95))

table(subset(res, exposure == "Lp(a)" & method == "Wald ratio")$pval < 0.05) %>% prop.test(p=c(0.95))


table(
	subset(res, method == "Weighted median")$pval < 0.05, 
	subset(res, method == "Weighted median")$exposure
)


table(
	subset(res, method == "MR Egger")$pval < 0.05, 
	subset(res, method == "MR Egger")$exposure
)



toggle_dev("test")
temp <- format_mr_results(res)
temp$fdr <- p.adjust(temp$pval, method="fdr")
temp$bonf <- p.adjust(temp$pval, method="bonferroni")



library(TwoSampleMR)

pbc_inst <- extract_instruments(287)
out <- extract_outcome_data(pbc_inst$SNP, 302)
dat <- harmonise_data(pbc_inst, out)

res <- mr(dat)


make_boxes <- function(dat, exposure_name, outcome_name)
{
	require(TwoSampleMR)
	require(gridExtra)
	require(ggplot2)
	require(ggrepel)

	temp <- subset(dat, outcome == outcome_name & exposure == exposure_name)
	temp$gene.exposure <- NA
	# temp$gene.exposure[! temp$gene.exposure %in% c("HMGCR", "PCSK9", "NPC1L1")] <- NA

	exposure_units <- temp$units.exposure[1]
	outcome_units <- temp$units.outcome[1]

	mrs <- mr_singlesnp(temp, all_method=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
	mrl <- mr_leaveoneout(temp)

	mrs$index <- 1:nrow(mrs)
	mrl$index <- 1:nrow(mrl)

	mrs <- merge(mrs, select(temp, SNP, gene.exposure), all.x=TRUE) %>% arrange(index)
	mrl <- merge(mrl, select(temp, SNP, gene.exposure), all.x=TRUE) %>% arrange(index)

	mrres <- mr(temp, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

	grid.arrange(
		mr_forest_plot(mrs)[[1]] + 
			labs(
				title="a)", 
				x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")")
			),
		mr_scatter_plot(mrres, temp)[[1]] + 
			labs(
				title="b)", 
				x=paste0("SNP effect on ", exposure_name),
				y=paste0("SNP effect on ", outcome_name)
			),
		mr_leaveoneout_plot(mrl)[[1]] + 
			labs(
				title="c)", 
				x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")"), 
				y="Excluded variant"
			),
		mr_funnel_plot(mrs)[[1]] + 
			labs(title="d)") +
			theme(legend.position="none"),
		ncol=2
	)
}

png(file="../images/sensitivity_boxes/pbc-triglycerides.png", width=900, height=900)
make_boxes(dat, dat$exposure[1], dat$outcome[1])
dev.off()



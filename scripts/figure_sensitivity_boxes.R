library(TwoSampleMR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(rmarkdown)

make_boxes <- function(dat, exposure_name, outcome_name)
{
	require(TwoSampleMR)
	require(ggplot2)
	require(ggrepel)

	temp <- subset(dat, outcome == outcome_name & exposure == exposure_name)
	temp$SNP <- paste0(temp$SNP, " (", temp$gene.exposure, ")")
	temp$gene.exposure[! temp$gene.exposure %in% c("HMGCR", "PCSK9", "NPC1L1")] <- NA

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
			) +
			geom_label_repel(aes(label=gene.exposure), point.padding = unit(0.5, "lines")),
		mr_leaveoneout_plot(mrl)[[1]] + 
			labs(
				title="c)", 
				x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")"), 
				y="Excluded variant"
			),
		mr_funnel_plot(mrs)[[1]] + 
			labs(title="d)") +
			theme(legend.position="none") +
			geom_label_repel(aes(label=gene.exposure), point.padding = unit(0.5, "lines")),
		ncol=2
	)
}

# Read in data
load("../results/lpa_ldl_trigs.RData")

# Get significant only and rename outcomes
outs <- select(outcomes, id, trait)
sig_res <- filter(res, method == "Inverse variance weighted", pval < 0.05, exposure %in% c("LDL cholesterol", "Triglycerides"))
sig_dat <- merge(subset(dat, outcome %in% sig_res$outcome), outs, by.x="id.outcome", by.y="id")
sig_dat$outcome <- sig_dat$trait
sig_res <- merge(sig_res, outs, by.x="id.outcome", by.y="id")
sig_res$outcome <- sig_res$trait
sig_res <- sig_res[order(sig_res$exposure, sig_res$outcome),]

# Make plots
dir.create("../images/sensitivity_boxes")

system("cp ../images/sensitivity_boxes/template.rmd ../images/sensitivity_boxes/supfig3.rmd")

i <- 1
for(i in 1:nrow(sig_res))
{
	fn1 <- paste0(gsub(" ", "_", sig_res$exposure[i]), "-", gsub(" ", "", gsub("[^[:alnum:] ]", "", sig_res$outcome[i])), ".png")
	fn <- paste0("../images/sensitivity_boxes/", fn1)
	print(fn)
	# png(file=fn, width=900, height=900)
	# make_boxes(sig_dat, outcome = sig_res$outcome[i], exposure = sig_res$exposure[i])
	# dev.off()
	nom <- paste0(sig_res$exposure[i], " on ", sig_res$outcome[i])
	system(paste0('echo "### ', nom, '" >> ../images/sensitivity_boxes/supfig3.rmd'))
	system(paste0('echo "" >> ../images/sensitivity_boxes/supfig3.rmd'))
	system(paste0('echo "![](', fn1, ')" >> ../images/sensitivity_boxes/supfig3.rmd'))
	system(paste0('echo "" >> ../images/sensitivity_boxes/supfig3.rmd'))

}

render("../images/sensitivity_boxes/supfig3.rmd")

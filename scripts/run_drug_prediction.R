library(TwoSampleMR)
library(ggplot2)
library(dplyr)

load("../results/drug_proxies.RData")

res <- filter(res, 
		outcome %in% c("Coronary heart disease || CARDIoGRAMplusC4D || 2015", "Type 2 diabetes || DIAGRAM || 2014", "Type 2 diabetes || DIAGRAMplusMetabochip || 2012"), 
		method == "Wald ratio"
	) %>%
	select(exposure, outcome, b, se) %>%
	mutate(outcome = sapply(strsplit(as.character(outcome), split="\\|"), function(x) { gsub(" $", "", x[[1]])}))

drug_trials <- read.csv("../data/drug_trials.csv") %>%
	mutate(b = b * 38.67459, se = se * 38.67459) %>%
	select(exposure, outcome, b, se)

res$test <- "Mendelian randomization"
drug_trials$test <- "Randomized controlled trial"
res <- rbind(res, drug_trials)

res$exposure <- factor(res$exposure, levels=c("HMGCR inhibitor", "NPC1L1 inhibitor", "PCSK9 inhibitor", "Lp(a) inhibitor", "APOC3 inhibitor"))
res$labels <- res$exposure
levels(res$labels) <- c(
	"Statins\nHMGCR variant\nTarget: LDL-C",
	"Ezetimibe\nNPC1L1 variant\nTarget: LDL-C",
	"Evolocumab\nPCSK9 variant\nTarget: LDL-C",
	"Apo(a) inhibitor\nLPA variant\nTarget: Lp(a)",
	"APOC3 inhibitor\nAPOC3 variant\nTarget: Triglycerides"
)


save(res, file="../results/drug_prediction.RData")

res$lor <- exp(res$b)
res$up <- exp(res$b + 1.96 * res$se)
res$lo <- exp(res$b - 1.96 * res$se)

subset(res, test=="Mendelian randomization", select=c(exposure, outcome, lor, up, lo))
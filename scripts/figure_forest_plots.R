library(TwoSampleMR)
library(ggplot2)
library(dplyr)
toggle_dev("test")

# Lipids on traits

load("../results/lpa_ldl_trigs.RData")
dir.create("../images/forest")

disease_plot <- subset(res, id.outcome %in% subset(outcomes, category == "Disease")$id) %>%
	forest_plot(., 
		by_category=TRUE,
		in_columns=TRUE,
		exponentiate=TRUE,
		trans="log2",
		xlab="Intervention on"
	)
ggsave(
	disease_plot, 
	file="../images/forest/disease_forest.png", 
	width=10, height=12
)

risk_factor_plot <- subset(res, id.outcome %in% subset(outcomes, category == "Risk factor" & !grepl("Intracranial volume", outcome))$id) %>%
	forest_plot(., 
		by_category=TRUE, 
		in_columns=TRUE,
		xlab="Intervention on",
		threshold = 0.05
	)
ggsave(
	risk_factor_plot, 
	file="../images/forest/risk_factor_forest.png", 
	width=10, height=16
)


# Drug targets on traits

load("../results/drug_proxies.RData")

disease_plot <- subset(res, id.outcome %in% subset(outcomes, category == "Disease")$id) %>%
	forest_plot(., 
		by_category=TRUE,
		in_columns=TRUE,
		exponentiate=TRUE,
		trans="log2",
		xlab="Predicted effect of\n",
		xlim=c(0.2, 3),
		threshold = 0.05
	)
ggsave(
	disease_plot, 
	file="../images/forest/drug_proxies_disease_forest.png", 
	width=14, height=9
)

risk_factor_plot <- subset(res, id.outcome %in% subset(outcomes, category == "Risk factor")$id) %>%
	forest_plot(., 
		by_category=TRUE, 
		in_columns=TRUE,
		xlab="Predicted effect of\n",
		xlim=c(-0.5, 0.5),
		threshold = 0.05
	)
ggsave(
	risk_factor_plot, 
	file="../images/forest/drug_proxies_risk_factor_forest.png", 
	width=14, height=14
)


library(ggplot2)

load("../results/drug_prediction.RData")

res$or <- exp(res$b)
res$or_lci <- exp(res$b - 1.96 * res$se)
res$or_uci <- exp(res$b + 1.96 * res$se)


ggplot(res, aes(y=test, x=or)) +
	geom_point(aes(colour=test), size=3) +
	scale_x_continuous(trans='log2',breaks=c(0.3,0.5,1,2,4)) + 
	geom_errorbarh(aes(xmin=or_lci, xmax=or_uci, colour=test), height=0) +
	facet_grid(labels ~ outcome, scale="free_x") +
	geom_vline(xintercept=1, linetype="dotted") + 
	theme(
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		strip.text.y=element_text(angle=0),
		legend.position="top",
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()) +
	labs(colour=NULL, y=NULL, x="Relative risk (95% CI) per SD decrease in lipid trait due to drug or genetic variant")
ggsave("../images/drug_prediction/mr_rct_comparison.png", width=8, height=4)


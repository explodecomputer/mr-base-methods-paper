library(TwoSampleMR)
library(ggrepel)
library(ggthemes)
toggle_dev("test")

## Volcano plot of lipids on traits

dir.create("../images/volcano/")

load("../results/lpa_ldl_trigs.RData")

fres <- format_mr_results(res)
fres$exposure[fres$exposure == "Lp(a)"] <- "Lipoprotein(a)"
fres <- subset(fres, ! outcome %in% c("Coronary heart disease", "Myocardial infarction") & !grepl("Intracranial volume", outcome))
fres$sig <- fres$pval < 0.05
fres$fdr <- p.adjust(fres$pval, method="fdr")

ggplot(fres, aes(x=effect, y=-log10(pval))) +
geom_vline(xintercept=0, linetype="dotted") +
geom_point(data=subset(fres, !sig)) +
geom_point(data=subset(fres, sig), aes(colour=category, size=fdr < 0.05)) +
facet_grid(. ~ exposure, scale="free") +
geom_label_repel(data=subset(fres, sig), aes(label=outcome, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=2, max.iter=3e3) +
# geom_text_repel(data=subset(fres, sig), aes(label=outcome, colour=category)) +
theme_bw() + 
theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      strip.text.x=element_text(size=16)
) +
# scale_colour_brewer(type="qual", palette="Dark2") +
# scale_fill_brewer(type="qual", palette="Dark2") +
labs(x="Effect on outcome (SD or log(OR)) per SD lowering of lipids", size="5% FDR") +
xlim(c(-0.7, 0.7))
ggsave("../images/volcano/volcano_boxes.png", width=16, height=8)



## Volcano plot of drug proxies on traits

load("../results/drug_proxies.RData")

fres <- format_mr_results(res)
fres <- subset(fres, ! outcome %in% c("Coronary heart disease", "Myocardial infarction") & category != "Lipid" & !grepl("Intracranial volume", outcome))
fres$sig <- fres$pval < 0.05

fres$fdr <- p.adjust(fres$pval, method="fdr")

ggplot(fres, aes(x=effect, y=-log10(pval))) +
geom_vline(xintercept=0, linetype="dotted") +
geom_point(data=subset(fres, !sig)) +
geom_point(data=subset(fres, sig), aes(colour=category, size=fdr < 0.05)) +
facet_grid(. ~ exposure, scale="free") +
geom_label_repel(data=subset(fres, sig), aes(label=outcome, fill=category), colour="white", segment.colour="black", point.padding = unit(0.7, "lines"), box.padding = unit(0.7, "lines"), segment.size=0.5, force=1, max.iter=3e3) +
theme_bw() + 
theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      strip.text.x=element_text(size=16)
) +
# scale_colour_brewer(type="qual", palette="Dark2") +
# scale_fill_brewer(type="qual", palette="Dark2") +
labs(x="Effect on outcome (SD or log(OR)) per SD lowering of lipids", size="5% FDR")
ggsave("../images/volcano/drug_proxies_volcano_boxes.png", width=19, height=8)

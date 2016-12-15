library(wordcloud)
set.seed(100)

load("../results/lpa_ldl_trigs.RData")

tab <- table(outcomes$subcategory)
colors <- colorRampPalette( brewer.pal( 8, "Dark2" ) )( length(tab) )
pdf("../images/schematic/subcategory_cloud.pdf")
wordcloud(names(tab), tab, scale=c(3,2), min.freq=1, random.order=T, rot.per=.150, colors=rev(colors))
dev.off()


temp <- data.frame(subcategory=unique(outcomes$subcategory), colors=colorRampPalette( brewer.pal( 8, "Dark2" ) )( length(unique(outcomes$subcategory)) ), stringsAsFactors=FALSE)
o <- merge(outcomes, temp, by="subcategory")

pdf("../images/schematic/trait_cloud.pdf")
wordcloud(o$trait,o$sample_size, scale=c(1.6,0.7), random.order=T, rot.per=.15, colors=o$colors)
dev.off()




library(dplyr)
library(GGally)
library(network)
library(sna)
library(ggplot2)
load("../results/lpa_ldl_trigs.RData")

# random graph
net = rgraph(10, mode = "graph", tprob = 0.5)
net = network(net, directed = FALSE)

net <- network(diag(57))

nom <- c(paste0("rs",1:57), "ldl", outcomes$trait)
net <- diag(length(nom))
net[1:57, 58] <- 1
net[58, 59:length(nom)] <- 1
net <- network(net, directed=TRUE)
network.vertex.names(net) = nom

net %v% "what" = c(rep("Instrument", 57), "Exposure", rep("Outcome", nrow(outcomes)))

ggnet2(net,arrow.size = 12, arrow.gap = 0.025, color = "what", mode="circle", palette="Set2")





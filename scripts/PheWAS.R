# install.packages("devtools")
# install_github("MRCIEU/TwoSampleMR")
# install.packages("xtable")
# install.packages("magrittr")
# install.packages("plyr")
# devtools::install_github("MRCIEU/MRInstruments")
library(devtools)
library(magrittr)
library(plyr)
library(TwoSampleMR)
library(MRInstruments)
data(gwas_catalog)

# SNP lookups for rs11065987, rs1250229 and rs4530754 in MR-Base and GWAS catalog. These were outlier SNPs identified in Mendelian randomization analyses of LDL cholesterol and coronary heart disease

########################
#SNP lookups in MR-Base#
########################

ao<-available_outcomes()
id<-ao$id
Dat <- extract_outcome_data(
    snps = c("rs11065987","rs1250229","rs4530754"), 
    outcomes = id
)

PheWAS1<-Dat 
PheWAS1<-PheWAS1[which(!duplicated(paste(PheWAS1$originalname.outcome,PheWAS1$SNP))),]

########################
#Lookups in GWAS catlog#
########################

PheWAS2<-gwas_catalog[gwas_catalog$SNP %in% c("rs11065987","rs1250229","rs4530754"),c("SNP","PubmedID","Phenotype","Phenotype_simple","beta","effect_allele","other_allele","units","Initial_sample_description","pval","eaf","date_added_to_MRBASE")]
PheWAS2<-PheWAS2[which(!duplicated(paste(PheWAS2$Phenotype_simple,PheWAS2$SNP))),]

###################################################################
# Format lookup results and merge GWAS catalog and MR-base lookups#
###################################################################

Ntests<-nrow(PheWAS2)+nrow(PheWAS1)
PheWAS1<-PheWAS1[which(PheWAS1$pval.outcome<0.05/Ntests),]
PheWAS2<-PheWAS2[which(PheWAS2$pval<0.05/Ntests),]
Excl1<-c("LDL cholesterol","Total cholesterol","HDL cholesterol","Myocardial infarction","Coronary heart disease","Ischemic stroke") #exclude vascular and lipid traits from lookups
PheWAS1$excl[!PheWAS1$originalname.outcome %in% Excl1]<-0
PheWAS1$excl[PheWAS1$originalname.outcome %in% Excl1]<-1
PheWAS1<-PheWAS1[order(PheWAS1$samplesize.outcome,decreasing=T),]
PheWAS1<-PheWAS1[!duplicated(PheWAS1$originalname.outcome),]

# Extract SD information for continuous traits
ao.sd<-ao[ao$id %in% PheWAS1$id.outcome,c("sd","id")]
PheWAS1<-merge(PheWAS1,ao.sd,by.x="id.outcome",by.y="id",all.x=T)
Pos<-grep("SD",PheWAS1$units.outcome,invert=T)
Pos2<-which(!is.na(PheWAS1$sd))
Pos<-Pos[Pos %in% Pos2]
PheWAS1$beta.outcome[Pos]<-PheWAS1$beta.outcome[Pos]/PheWAS1$sd[Pos]
PheWAS1$se.outcome[Pos]<-PheWAS1$se.outcome[Pos]/PheWAS1$sd[Pos]

PheWAS1$source<-"MRBASE"
Excl2<-c("Cholesterol, total (unit increase)","Cholesterol, total (mg/dL decrease)","Cholesterol, total (unit decrease)","Coronary artery disease or ischemic stroke","Coronary artery disease or large artery stroke","LDL cholesterol (mg/dL decrease)","LDL cholesterol (unit decrease)")
PheWAS2$excl[!PheWAS2$Phenotype %in% Excl2]<-0
PheWAS2$excl[PheWAS2$Phenotype %in% Excl2]<-1
PheWAS2$source<-"GWAScatalog"
names(PheWAS1)[names(PheWAS1) %in% c("beta.outcome","se.outcome","samplesize.outcome","pval.outcome","eaf.outcome","effect_allele.outcome","other_allele.outcome","units.outcome","consortium.outcome","pmid.outcome","originalname.outcome","proxy.outcome")]<-c("beta","se","samplesize","pval","eaf","effect_allele","other_allele","units","consortium","pmid","exposure","proxy")

names(PheWAS2)[names(PheWAS2) %in% c("PubmedID","Phenotype_simple","Initial_sample_description")] <- c("pmid","exposure","samplesize")


PheWAS<-rbind.fill(PheWAS1,PheWAS2)
PheWAS<-PheWAS[,c("SNP","beta","se","samplesize","pval","eaf","effect_allele","other_allele","units","sd","consortium","pmid","exposure","proxy","date_added_to_MRBASE","source","excl")]

# Flip all SNPs to reflect the LDL increasing allele
Pos<-which(PheWAS$SNP=="rs11065987" & PheWAS$effect_allele =="G")

PheWAS$beta[Pos]<-PheWAS$beta[Pos]*-1
PheWAS$eaf[Pos]<-1-PheWAS$eaf[Pos]
other_allele<-PheWAS$effect_allele[Pos]
effect_allele<-PheWAS$other_allele[Pos]
PheWAS$other_allele[Pos]<-other_allele
PheWAS$effect_allele[Pos]<-effect_allele

Pos<-which(PheWAS$SNP=="rs1250229" & PheWAS$effect_allele =="T")
PheWAS$beta[Pos]<-PheWAS$beta[Pos]*-1
PheWAS$eaf[Pos]<-1-PheWAS$eaf[Pos]
other_allele<-PheWAS$effect_allele[Pos]
effect_allele<-PheWAS$other_allele[Pos]
PheWAS$other_allele[Pos]<-other_allele
PheWAS$effect_allele[Pos]<-effect_allele

Pos<-which(PheWAS$SNP=="rs4530754" & PheWAS$effect_allele =="G")
PheWAS$beta[Pos]<-PheWAS$beta[Pos]*-1
PheWAS$eaf[Pos]<-1-PheWAS$eaf[Pos]
other_allele<-PheWAS$effect_allele[Pos]
effect_allele<-PheWAS$other_allele[Pos]
PheWAS$other_allele[Pos]<-other_allele
PheWAS$effect_allele[Pos]<-effect_allele
PheWAS$units<-gsub("increase","",PheWAS$units)
PheWAS$units<-gsub("decrease","",PheWAS$units)

PheWAS<-PheWAS[order(PheWAS$excl,PheWAS$exposure),]

write.table(PheWAS,"snplookups.txt",sep="\t",col.names=T,row.names=F,quote=F)

#####################################################################################################################################
# Further Mendelian randomization analyses were conducted, to assess effect of the traits identified above on coronary heart disease
#####################################################################################################################################
# First, we must define instruments using either MR-Base or the GWAS catalog 

id<-PheWAS1$id.outcome[grep("pressure",PheWAS1$exposure,invert=T)] 
exposure_dat1 <- extract_instruments(id)

# Define instruments using GWAS catalog for traits retrieved from lookups and which can't be defined using MR-Base:
# The following traits had insufficient data for instrument extraction in MR-Base: 
	# c("Diastolic blood pressure" "Systolic blood pressure" , "Mean arterial pressure","Hematocrit (% decrease)","Tetralogy of Fallot") 

# hematocrit
ht_gwas<-gwas_catalog[which(gwas_catalog$PubmedID == 19862010),]
ht_gwas$Phenotype<-ht_gwas$Phenotype_simple 
ht_gwas<-ht_gwas[ht_gwas$Phenotype %in% c("Hematocrit"),]
ht_gwas$units[ht_gwas$Phenotype== "Hematocrit"]<-"%"

# systolic and diastolic blood pressure
bp_gwas1<-gwas_catalog[which(gwas_catalog$PubmedID == 21909115),]
bp_gwas1$Phenotype<-bp_gwas1$Phenotype_simple 
bp_gwas1<-bp_gwas1[bp_gwas1$Phenotype_simple != "Hypertension",]
bp_gwas1$units<-"mmHg"

# mean arterial pressure
bp_gwas2<-gwas_catalog[which(gwas_catalog$PubmedID == 21909110),]
bp_gwas2<-bp_gwas2[bp_gwas2$MAPPED_TRAIT_EFO =="mean arterial pressure",]
bp_gwas2$Phenotype<-bp_gwas2$MAPPED_TRAIT_EFO
bp_gwas2$units<-"mmHg"
bp_gwas<-rbind(bp_gwas1,bp_gwas2)

# Tetralogy of Fallot
tet_gwas1<-gwas_catalog[which(gwas_catalog$PubmedID == 23297363),]
tet_gwas1$units<-"log odds"

gwas<-do.call(rbind,list(tet_gwas1,bp_gwas,ht_gwas))
exposure_dat2<-format_data(gwas)
exposure_dat2$data_source.exposure<-"GWAS catalog"

exposure_dat<-rbind.fill(exposure_dat1,exposure_dat2)
outcome_dat <- extract_outcome_data(unique(exposure_dat$SNP), 7, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results_dat <- mr(dat)
mr_results<-mr_results_dat
# extract SD information and merge with MR results
ao.sd<-ao[ao$id %in% mr_results$id.exposure,c("sd","id")]
mr_results<-merge(mr_results,ao.sd,by.x="id.exposure",by.y="id",all.x=T)
Pos<-grep("SD",mr_results$exposure,invert=T)
Pos2<-which(!is.na(mr_results$sd))
Pos<-Pos[Pos %in% Pos2]

#run Mendelian randomization analyses
mr_results$b[Pos]<-mr_results$b[Pos]*mr_results$sd[Pos]
mr_results$se[Pos]<-mr_results$se[Pos]*mr_results$sd[Pos]
mr_results$lci<-mr_results$b-1.96*mr_results$se
mr_results$uci<-mr_results$b+1.96*mr_results$se
mr_ivw<-subset_on_method(mr_results)
mr_ivw$fdr <- p.adjust(mr_ivw$pval, method="fdr")
mr_ivw$bonf<-5e-2/nrow(mr_ivw)
write.table(mr_results,"mr_results.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(mr_ivw,"mr_ivw.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Infer missing standard errors from Z score and effect size
# Hemotocrit SE missing for rs11065987
p<-1.00E-12
b<-0.17
z<-qnorm(p/2,lower.tail=F)
se<-b/z

# Tetrology of fallot for rs11065987
b<--0.296 
p<-8.00E-11
z<-qnorm(p/2,lower.tail=F)
se<-abs(b)/z



# install.packages("devtools")
# library(devtools)
# install_github("MRCIEU/TwoSampleMR")
# install.packages("plyr")
setwd("~/mr-base-methods-paper/data/")
library(plyr)
library(TwoSampleMR)

# SNP lookups for rs11065987, rs1250229 and rs4530754 in MR-Base and GWAS catalog. These were outlier SNPs identified in Mendelian randomization analyses of LDL cholesterol and coronary heart disease

########################
#SNP lookups in MR-Base#
########################

load("gwas_catalog_Nov2017.RData")
load("available_outcomes_Nov2017.RData")

id<-ao$id
Dat <- extract_outcome_data(
    snps = c("rs11065987","rs1250229","rs4530754"), 
    outcomes = id
)

PheWAS1<-Dat[order(Dat$samplesize.outcome,decreasing=T),]
PheWAS1<-PheWAS1[which(!duplicated(paste(PheWAS1$originalname.outcome,PheWAS1$SNP))),]

########################
#Lookups in GWAS catlog#
########################

PheWAS2<-gwas_catalog[gwas_catalog$SNP %in% c("rs11065987","rs1250229","rs4530754"),c("SNP","PubmedID","Phenotype","Phenotype_simple","beta","effect_allele","other_allele","units","Initial_sample_description","pval","eaf","date_added_to_MRBASE")]

Sample.size<-regmatches(PheWAS2$Initial_sample_description,gregexpr("[[:digit:]]+\\,*[[:digit:]]*",PheWAS2$Initial_sample_description))

PheWAS2$samplesize<-unlist(lapply(1:length(Sample.size),FUN=function(x) sum(as.numeric(gsub(",","",unlist(Sample.size[x]))))))
PheWAS2<-PheWAS2[order(PheWAS2$samplesize,decreasing=T),]
PheWAS2<-PheWAS2[which(!duplicated(paste(PheWAS2$Phenotype_simple,PheWAS2$SNP))),]

	
###################################################################
# Format lookup results and merge GWAS catalog and MR-base lookups#
###################################################################

Ntests<-nrow(PheWAS2)+nrow(PheWAS1)
PheWAS1<-PheWAS1[which(PheWAS1$pval.outcome<0.05/Ntests),]
PheWAS2<-PheWAS2[which(PheWAS2$pval<0.05/Ntests),]
Excl1<-c("LDL cholesterol","Total cholesterol","HDL cholesterol","Myocardial infarction","Coronary heart disease","Ischemic stroke") #exclude vascular disease and lipid traits from lookups
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

PheWAS1$phewas.source<-"MRBASE"
Excl2<-c("Cholesterol, total (unit increase)","Cholesterol, total (mg/dL decrease)","Cholesterol, total (unit decrease)","Coronary artery disease or ischemic stroke","Coronary artery disease or large artery stroke","LDL cholesterol (mg/dL decrease)","LDL cholesterol (unit decrease)")
PheWAS2$excl[!PheWAS2$Phenotype %in% Excl2]<-0
PheWAS2$excl[PheWAS2$Phenotype %in% Excl2]<-1
PheWAS2$phewas.source<-"GWAScatalog"
names(PheWAS1)[names(PheWAS1) %in% c("beta.outcome","se.outcome","samplesize.outcome","pval.outcome","eaf.outcome","effect_allele.outcome","other_allele.outcome","units.outcome","consortium.outcome","pmid.outcome","originalname.outcome","proxy.outcome")]<-c("beta","se","samplesize","pval","eaf","effect_allele","other_allele","units","consortium","pmid","exposure","proxy")

names(PheWAS2)[names(PheWAS2) %in% c("PubmedID","Phenotype_simple")] <- c("pmid","exposure")


PheWAS<-rbind.fill(PheWAS1,PheWAS2)
PheWAS<-PheWAS[,c("SNP","beta","se","samplesize","pval","eaf","effect_allele","other_allele","units","sd","consortium","pmid","exposure","proxy","date_added_to_MRBASE","phewas.source","excl")]

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

#####################################################################################################################################
# Further Mendelian randomization analyses were conducted, to assess effect of the traits identified above on coronary heart disease
#####################################################################################################################################
# First, we must define instruments using either MR-Base or the GWAS catalog 

id<-PheWAS1$id.outcome 
exposure_dat1 <- extract_instruments(id)
exposure_dat1<-exposure_dat1[exposure_dat1$SNP!="rs11065987",]
exposure_dat1$instrument.source<-"MR_Base"

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
# bp_gwas2$sd<-12.8
bp_gwas<-rbind(bp_gwas1,bp_gwas2)

# SD for mean arterial pressure = 12.8

# Tetralogy of Fallot
tet_gwas1<-gwas_catalog[which(gwas_catalog$PubmedID == 23297363),]
tet_gwas1$units<-"log odds"

# Primary biliary cirrhosis
pbc_gwas<-gwas_catalog[which(gwas_catalog$PubmedID == 26394269),]
pbc_gwas$units<-"log odds"

gwas<-do.call(rbind.fill,list(tet_gwas1,bp_gwas,ht_gwas,pbc_gwas))
exposure_dat2<-format_data(gwas)
exposure_dat2<-exposure_dat2[exposure_dat2$SNP!="rs11065987",]
exposure_dat3<-clump_data(exposure_dat2)
exposure_dat2$instrument.source<-"GWAS catalog"

exposure_dat<-rbind.fill(exposure_dat1,exposure_dat2)

outcome_dat <- extract_outcome_data(unique(exposure_dat$SNP), 7, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results_dat <- mr(dat,method_list = c("mr_wald_ratio","mr_ivw_mre","mr_ivw_fe"))
mr_results<-mr_results_dat
mr_results<-merge(mr_results,exposure_dat[,c("id.exposure","instrument.source")],by="id.exposure")
# extract SD information and merge with MR results
ao.sd<-ao[ao$id %in% mr_results$id.exposure,c("sd","id","unit")]
mr_results<-merge(mr_results,ao.sd,by.x="id.exposure",by.y="id",all.x=T)
mr_results$sd[mr_results$exposure=="Hematocrit"]<-3.339
mr_results$unit[mr_results$exposure=="Hematocrit"]<-"%"
mr_results$sd[mr_results$exposure=="Systolic blood pressure"]<-18.2
mr_results$unit[mr_results$exposure=="Systolic blood pressure"]<-"mmHg"
mr_results$sd[mr_results$exposure=="Diastolic blood pressure"]<-10.7
mr_results$unit[mr_results$exposure=="Diastolic blood pressure"]<-"mmHg"
mr_results$sd[mr_results$exposure=="mean arterial pressure"]<-12.8
mr_results$unit[mr_results$exposure=="mean arterial pressure"]<-"mmHg"
Pos<-grep("SD",mr_results$unit,invert=T)
Pos2<-which(!is.na(mr_results$sd))
Pos<-Pos[Pos %in% Pos2]

mr_results$b[Pos]<-mr_results$b[Pos]*mr_results$sd[Pos]
mr_results$se[Pos]<-mr_results$se[Pos]*mr_results$sd[Pos]
mr_results<-mr_results[order(mr_results$se,decreasing=T),]
mr_ivw<-mr_results[!duplicated(mr_results$exposure),] #drop duplicate results, keeping the result that had the larger standard error. The random effects model usually gives a larger standard error except when there is underdisperson in the SNP effects on the outcome, in which case a fixed effects model is perferable for the standard error



# create suppl table 4

mr_ivw<-mr_ivw[!mr_ivw$exposure %in% c("Ischemic stroke || id:1108","Myocardial infarction || id:798","Total cholesterol || id:301","HDL cholesterol || id:299"),] #excluding lipid and vascular disease traits other than LDL cholesterol and coronary heart disease

PheWAS<-PheWAS[PheWAS$SNP=="rs11065987",]
PheWAS<-PheWAS[!PheWAS$exposure %in% c("Total cholesterol","HDL cholesterol","Myocardial infarction","Ischemic stroke", "Cholesterol, total","Coronary artery disease or ischemic stroke","Coronary artery disease or large artery stroke"),]
PheWAS$exposure<-gsub("Hemoglobin","Haemoglobin concentration",PheWAS$exposure)

PheWAS<-PheWAS[order(PheWAS$phewas.source,decreasing=T),]
PheWAS<-PheWAS[!duplicated(PheWAS$exposure),] #get rid of duplicates for BMI and LDL cholesterol, owing to lookups in GWAS catalog and MR-Base. Only lookups in MR-Base had non-missing standard errors
names(PheWAS)[names(PheWAS)=="beta"]<-"b"
names(PheWAS)[names(PheWAS)=="units"]<-"unit"
mr_ivw<-split_outcome(mr_ivw)
mr_ivw<-split_exposure(mr_ivw)

mr_ivw$exposure<-gsub("mean arterial pressure","Mean arterial pressure",mr_ivw$exposure)
mr_ivw_phewas<-merge(mr_ivw,PheWAS[,c("exposure","b","se","pval","unit","pmid","phewas.source","effect_allele","other_allele")],by="exposure")
mr_ivw_phewas$x<-"effect of trait on CHD"
mr_ivw_phewas$y<-"effect of SNP on trait"
mr_ivw_phewas$z<-"effect of trait on CHD scaled to magnitude of SNP-trait effect"
mr_ivw_phewas$b.y[mr_ivw_phewas$phewas.source=="GWAScatalog" & !is.na(mr_ivw_phewas$sd)]<-mr_ivw_phewas$b.y[mr_ivw_phewas$phewas.source=="GWAScatalog" & !is.na(mr_ivw_phewas$sd)]/mr_ivw_phewas$sd[mr_ivw_phewas$phewas.source=="GWAScatalog" & !is.na(mr_ivw_phewas$sd)]
mr_ivw_phewas$b.z <- mr_ivw_phewas$b.x * mr_ivw_phewas$b.y #scale the effect of the trait on CHD by the magnitude of the SNP-trait effect
mr_ivw_phewas$se.z <- abs(mr_ivw_phewas$se.x * mr_ivw_phewas$b.y) #scale the effect of the trait on CHD by the magnitude of the SNP-trait effect
mr_ivw_phewas$b.x[mr_ivw_phewas$exposure=="Coronary heart disease"]<-NA
mr_ivw_phewas$se.x[mr_ivw_phewas$exposure=="Coronary heart disease"]<-NA
mr_ivw_phewas$pval.x[mr_ivw_phewas$exposure=="Coronary heart disease"]<-NA
mr_ivw_phewas$b.z[mr_ivw_phewas$exposure=="Coronary heart disease"]<-NA
mr_ivw_phewas$se.z[mr_ivw_phewas$exposure=="Coronary heart disease"]<-NA
mr_ivw_phewas$pval.z[mr_ivw_phewas$exposure=="Coronary heart disease"]<-NA

trait.order<-c("Coronary heart disease","LDL cholesterol","Body mass index","Hip circumference","Systolic blood pressure","Diastolic blood pressure","Mean arterial pressure","Red blood cell count","Haemoglobin concentration","Hematocrit","Packed cell volume","Platelet count","C-glycosyltryptophan*","Erythronate*","Kynurenine","Serum cystatin C (eGFRcys)","Urate","Inflammatory bowel disease","Crohn's disease","Rheumatoid arthritis","Primary biliary cirrhosis","Celiac disease","Tetralogy of Fallot")
mr_ivw_phewas<-mr_ivw_phewas[match(trait.order,mr_ivw_phewas$exposure),]
mr_ivw_phewas$exposure<-gsub("*","",mr_ivw_phewas$exposure)

# Infer missing standard errors for SNP-trait effect (se.y) from P value and effect size (pval.y and b.y)
# Hemotocrit SE missing for rs11065987
Pos<-is.na(mr_ivw_phewas$se.y)
p<-mr_ivw_phewas$pval.y[Pos] 
# 1.00E-12
b<-mr_ivw_phewas$b.y[Pos] 
z<-qnorm(as.numeric(p)/2,lower.tail=F)
mr_ivw_phewas$se.y[Pos]<-abs(b/z)
#round off betas, standard errors and pvalues
mr_ivw_phewas$b.x<-round(mr_ivw_phewas$b.x,3)
mr_ivw_phewas$b.y<-round(mr_ivw_phewas$b.y,3)
mr_ivw_phewas$b.z<-round(mr_ivw_phewas$b.z,3)
mr_ivw_phewas$se.x<-round(mr_ivw_phewas$se.x,3)
mr_ivw_phewas$se.y<-round(mr_ivw_phewas$se.y,3)
mr_ivw_phewas$se.z<-round(mr_ivw_phewas$se.z,3)
mr_ivw_phewas$pval.x<-formatC(mr_ivw_phewas$pval.x, format = "e", digits = 2)
mr_ivw_phewas$pval.y<-formatC(mr_ivw_phewas$pval.y, format = "e", digits = 2)
mr_ivw_phewas$method.note<-NA
mr_ivw_phewas$method.note[mr_ivw_phewas$method=="Inverse variance weighted (fixed effects)"]<-"se.x estimated using fixed effects model because of underdispersion in causal estimates amongst SNPs"

mr_ivw_phewas[is.na(mr_ivw_phewas$se.y),c("exposure","b.y","se.y","pval.y","sd")]

names(mr_ivw_phewas)[names(mr_ivw_phewas)=="nsnp"]<-"nsnp.instrument"

mr_ivw_phewas$b.se.z<-paste(mr_ivw_phewas$b.z," (",mr_ivw_phewas$se.z,")",sep="")
mr_ivw_phewas$b.se.y<-paste(mr_ivw_phewas$b.y," (",mr_ivw_phewas$se.y,")",sep="")
write.table(mr_ivw_phewas[,c("exposure","b.se.y","b.se.z","nsnp.instrument","b.y","se.y","pval.y","pmid","b.x","se.x","b.z","se.z","pval.x","phewas.source","instrument.source","method","method.note","x","y","z")],"~/mr_base_paper_elifev2/results/mr_ivw_phewas_clump2.txt",sep="\t",col.names=T,row.names=F,quote=F)

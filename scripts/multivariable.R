library(TwoSampleMR)
toggle_dev("test")
load("../results/lpa_ldl_trigs.RData")

pbc_multi <- multivariable_mr(c(299,300,302), 287)
t2d_multi <- multivariable_mr(c(299,300,302), 23)
crea_multi <- multivariable_mr(c(299,300,302), 1104)

l <- list()
for(i in 1:nrow(outcomes))
{
	message(i, "of", nrow(outcomes))
	l[[i]] <- multivariable_mr(c(299,300,302), outcomes$id[i])
}


# data <- read_excel("~/Google Drive/working space/Postdoc_year1_year2/Lipids-BMD-MR/SNPsnap_lipids_185/lipids-instruments-Do-185.xlsx",1)
# data[159,9]<-0
# data[185,13]<-0
# data<-data[1:184,]

# #read CHD beta and se 
# by<-data$CAD_beta
# sigmay<-data$CAD_se

# #or read bmd beta and se 
# by<-data$BETA_BMD
# sigmay<-data$SE_BMD

# #method 1 - Do's version (Gib is using)
# beta1 = lm(lm(by ~ bx2 + bx3)$res ~ bx1)$coef[2]
# beta2 = lm(lm(by ~ bx1 + bx3)$res ~ bx2)$coef[2]
# beta3 = lm(lm(by ~ bx1 + bx2)$res ~ bx3)$coef[2]
# betas <- c(beta1,beta2,beta3)
# beta1se = summary(lm(lm(by ~ bx2 + bx3)$res ~ bx1))$coef[2,2]
# beta2se = summary(lm(lm(by ~ bx1 + bx3)$res ~ bx2))$coef[2,2]
# beta3se = summary(lm(lm(by ~ bx1 + bx2)$res ~ bx3))$coef[2,2]
# Std.errs <-c(beta1se,beta2se,beta3se)
# zscore <- betas / Std.errs
# pval <- 2*pnorm(-abs(zscore))
# results.old.version <- cbind(betas,Std.errs,zscore,pval)
# row.names(results.old.version)<-c("LDL-c","HDL-c","TG")

# #method 2 - constrain the intercept - Stephen's version
# betas <- summary(lm(by~bx1+bx2+bx3-1, weights=sigmay^-2))$coef[,1]
# Std.errs <- summary(lm(by~bx1+bx2+bx3-1, weights=sigmay^-2))$coef[,2]/summary(lm(by~bx1+bx2+bx3-1, weights=sigmay^-2))$sigma
# zscore <- betas / Std.errs
# pval <- 2*pnorm(-abs(zscore))
# results.no.intercept <- cbind(betas,Std.errs,zscore,pval)


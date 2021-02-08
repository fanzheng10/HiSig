library(pROC)

setwd("~/Desktop/UCSD_Research/2021/hisig_paper/N3_m5_p0.3_1")
results <- read.delim('gsea_results.txt')
pos = readLines('selected_terms.txt')
results$sel = 0
results[results$V1 %in% pos,]$sel= 1
roc(results$sel ~ results$p,plot=TRUE,print.auc=TRUE,col="green",print.auc.y=0.3,
    lwd =2,legacy.axes=TRUE,main="ROC Curves")

hisig <- read.delim('sim_ms_impact_summary.tsv')
results$hisig_p = 1
hisig <- hisig[order(hisig$System_name),]
results[results$V1 %in% hisig$System_name,]$hisig_p = hisig$p
roc(results$sel ~ results$hisig_p,plot=TRUE,print.auc=TRUE,col="red",print.auc.y=0.5,
    lwd =2,legacy.axes=TRUE,add=T)

results$lasso = 0
lasso = read.delim('normal_lasso_no_identity.txt')
results[results$V1 %in% lasso$V1,]$lasso = lasso$coef

roc(results$sel ~ results$lasso,plot=TRUE,print.auc=TRUE,col="blue",print.auc.y=0.4,
    lwd =2,legacy.axes=TRUE,add=T)

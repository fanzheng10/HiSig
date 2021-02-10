library(ggplot2)
library(ggpubr)

setwd("~/Desktop/UCSD_Research/2021/hisig_paper/N3_m5_p0.3_1")

coef = read.table('sim_ms_impact.coef', header = F, row.names = 1)
lambdas = coef[1,]
coef = as.matrix(coef[2:nrow(coef),])
coef1 = coef[,which(colSums(coef)>0)]
lambdas = as.numeric(lambdas[,which(colSums(coef)>0)])
coef = coef1
coef = scale(coef, center=FALSE, scale=colSums(coef))

gene_features = as.numeric(colSums(coef[1:625,]))
L3_features = as.numeric(colSums(coef[626:630,]))
L2_features = as.numeric(colSums(coef[631:655,]))
L1_features = as.numeric(colSums(coef[656:780,]))

df_gene = data.frame(features=gene_features, lambdas)
df_gene$layer = 'gene'
df_L3 = data.frame(features=L3_features, lambdas)
df_L3$layer = 'L3'
df_L2 = data.frame(features=L2_features, lambdas)
df_L2$layer = 'L2'
df_L1 = data.frame(features=L1_features, lambdas)
df_L1$layer = 'L1'

df_all = rbind(df_gene, df_L3, df_L2, df_L1)

type_order = c('gene', 'L1', 'L2', 'L3')
df_all$layer = factor(df_all$layer, levels= type_order)
# df_all$lambdas = factor(df_all$lambdas)
p<-ggplot(df_all, aes(fill=layer, y=features, x=lambdas)) +
  geom_bar(position="fill", stat="identity",
           width=0.02,
           alpha=1,
  ) +theme_classic(base_size = 20) + scale_x_log10() +
  scale_fill_manual(values=rev(c('grey25', 'grey50', 'grey75', 'grey90'))) +
  # theme(axis.ticks=element_text(size=14),
  #       axis.title=element_text(size=18)) +
  xlab('Lambda') + ylab('Relative contribution\nof feature')

ggsave('../HiSig_vignette/simulation/fig/relative_contribution.pdf', p, device='pdf', width=7.5, height=5)


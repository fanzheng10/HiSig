library(viper)

data(bcellViper, package="bcellViper")

signature <- rowTtest(dset, "description", c("CB", "CC"), "N")
signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
write.table(signature, 'signature.txt',row.names = F, col.names = F, quote = F)
nullmodel <- ttestNull(dset, "description", c("CB", "CC"), "N", per = 1000, repos = TRUE, verbose = FALSE)
#
mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE)

# write regulon

# tf = names(regulon)
# write.table(tf, 'tf.txt', row.names = F, col.names = F, quote = F)
#
# # write targets
#
# targets = c()
# for (i in 1:length(regulon)) {
#   targets = c(targets, names(regulon[[i]]$tfmode))
# }
# targets = sort(unique(targets))
#
# write.table(targets, 'target.txt', row.names=F, col.names=F, quote=F)
#
# # write interactions
# tf_ids = c()
# tg_ids = c()
#
# for (i in 1:length(regulon)) {
#   tgs = names(regulon[[i]]$tfmode)
#   tf_id = rep.int(i, length(tgs))
#   tg_id = match(tgs, targets)
#   tf_ids = c(tf_ids, tf_id)
#   tg_ids = c(tg_ids, tg_id)
# }
#
# conn = cbind(tg_ids, tf_ids)
# write.table(conn, 'conn.txt', row.names = F, col.names=F, quote=F) # note this starts from 1

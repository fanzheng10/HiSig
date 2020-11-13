library(HiSig)

genes = readLines('target.txt')
terms = readLines('tf.txt')
data <- load_data('conn.txt', 'response.txt', genes=genes, terms=terms)

hisig_out <- hisig_fit(data)
beta_max_rand <- hisig_fit_rand(data, hisig_out$lambda, batch=100)

impact <- cbind(hisig_out$beta.max, beta_max_rand)

df <- parse_hisig(data, impact, term.names = terms, gene.names = genes)
output(df, 'hisig_summary.txt')

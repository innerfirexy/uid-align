# analyze the swbd_dysf_adjacent_pairs_count.txt results
# Yang Xu
# 10/3/2016

library(data.table)
library(ggplot2)

df = read.table(file = 'swbd_dysf_adjacent_pairs_count.txt', sep = ',', header = F)
dt = data.table(df)

colnames(dt) = c('rule', 'freq')
setkey(dt, freq)
setorder(dt, -freq)

# plot
plot(density(dt$freq))

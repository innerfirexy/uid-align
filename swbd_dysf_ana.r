# analyze the swbd_dysf_adjacent_pairs_count.txt results
# Yang Xu
# 10/3/2016

library(data.table)
library(ggplot2)
# library(reshape)

df = read.table(file = 'swbd_dysf_adjacent_pairs_count.txt', sep = ',', header = F)
dt = data.table(df)

colnames(dt) = c('rule', 'repeatFreq')
setkey(dt, repeatFreq)
setorder(dt, -repeatFreq)

# plot
plot(density(dt$repeatFreq))


# read data of all subrules freq
df2 = read.table(file = 'all_subrules_freq.txt', sep = ',', header = F)
dt2 = data.table(df2)
colnames(dt2) = c('rule', 'priorFreq')
setkey(dt2, priorFreq)
setorder(dt2, -priorFreq)

# join dt and dt2
setkey(dt, rule)
setkey(dt2, rule)
dt.join = dt[dt2, nomatch = 0]
setorder(dt.join, -repeatFreq)

# plot dt.join
dt.join.sample = dt.join[1:10,]
dt.join.sample$rule = as.character(dt.join.sample$rule)
dt.join.sample$rule = factor(dt.join.sample$rule, levels = dt.join.sample$rule)

dt.join.p = melt(dt.join.sample, id.vars = 'rule')

p1 = ggplot(dt.join.p, aes(x = rule, y = value, group = variable)) +
    geom_point(aes(shape = variable))


## next, compute probBoost for each rule

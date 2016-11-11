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
colnames(dt2) = c('rule', 'priorCount')
setkey(dt2, priorCount)
setorder(dt2, -priorCount)
dt2$priorFreq = dt2$priorCount / sum(dt2$priorCount)

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
df3 = read.table(file = 'pb_results_10.txt', sep = ',', header = F)
colnames(df3) = c('rule', 'probBoost')

dt.ss = dt.join.sample[, .(rule, priorFreq)]
dt.ss$probBoost = df3$probBoost
dt.ss.p = melt(dt.ss, id.vars = 'rule')

p = ggplot(dt.ss.p, aes(x = rule, y = value, group = variable)) +
    geom_point(aes(shape = variable, color = variable)) +
    geom_line(aes(lty = variable, color = variable)) +
    scale_linetype_manual(values = c(2, 1),
        name = 'Type', breaks = c('priorFreq', 'probBoost'),
        labels = c('prior probability', 'probability boost')) +
    scale_color_discrete(name = 'Type', breaks = c('priorFreq', 'probBoost'),
        labels = c('prior probability', 'probability boost')) +
    scale_shape_discrete(name = 'Type', breaks = c('priorFreq', 'probBoost'),
        labels = c('prior probability', 'probability boost')) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position=c(.4, .8))

pdf('probBoost_vs_prior_10.pdf', 5, 5)
plot(p)
dev.off()

# include the most frequent 20 rules
setorder(dt2, -priorCount)
dt.ss2 = dt2[1:20, .(rule, priorFreq)]
dt.ss2$rule = as.character(dt.ss2$rule)
dt.ss2$rule = factor(dt.ss2$rule, levels = dt.ss2$rule)

df4 = read.table(file = 'pb_results_20.txt', sep = ',', header = F)
colnames(df4) = c('rule', 'probBoost')
dt.ss2$probBoost = df4$probBoost
dt.ss2.p = melt(dt.ss2, id.vars = 'rule')

p = ggplot(dt.ss2.p, aes(x = rule, y = value, group = variable)) +
    geom_point(aes(shape = variable, color = variable)) +
    geom_line(aes(lty = variable, color = variable)) +
    scale_linetype_manual(values = c(2, 1),
        name = 'Type', breaks = c('priorFreq', 'probBoost'),
        labels = c('prior probability', 'probability boost')) +
    scale_color_discrete(name = 'Type', breaks = c('priorFreq', 'probBoost'),
        labels = c('prior probability', 'probability boost')) +
    scale_shape_discrete(name = 'Type', breaks = c('priorFreq', 'probBoost'),
        labels = c('prior probability', 'probability boost')) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position=c(.2, .8))

pdf('probBoost_vs_prior_20.pdf', 10, 5)
plot(p)
dev.off()

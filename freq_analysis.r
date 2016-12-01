# Analyze the preliminary frequency data in resutls folder
# Yang Xu
# 11/30/2016

library(data.table)
library(ggplot2)

# load data
dt1 = fread('results/NP -> DT NN_freq.txt')
dt2 = fread('results/NP -> NN_freq.txt')
dt3 = fread('results/NP -> NNS_freq.txt')

setnames(dt1, c('prime', 'target'))
setnames(dt2, c('prime', 'target'))
setnames(dt3, c('prime', 'target'))

dt1.plot = melt(dt1)
dt2.plot = melt(dt2)
dt3.plot = melt(dt3)

# plot
p1 = ggplot(dt1.plot, aes(x = variable, y = value)) +
    stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', width = .2) +
    stat_summary(fun.y = mean, geom = 'point') +
    labs(x = 'word position', y = 'frequency') + ggtitle('NP=DT+NN')
pdf('results/NP=DT+NN_freq.pdf', 5, 5)
plot(p1)
dev.off()


p2 = ggplot(dt2.plot, aes(x = variable, y = value)) +
    stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', width = .2) +
    stat_summary(fun.y = mean, geom = 'point') +
    labs(x = 'word position', y = 'frequency') + ggtitle('NP=NN')
pdf('results/NP=NN_freq.pdf', 5, 5)
plot(p2)
dev.off()

p3 = ggplot(dt3.plot, aes(x = variable, y = value)) +
    stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', width = .2) +
    stat_summary(fun.y = mean, geom = 'point') +
    labs(x = 'word position', y = 'frequency') + ggtitle('NP=NNS')
pdf('results/NP=NNS_freq.pdf', 5, 5)
plot(p3)
dev.off()

t.test(dt$prime, dt$target) # n.s.

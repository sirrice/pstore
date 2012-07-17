## #
## #
## # Compare write costs (overhead) with query costs
## #
## #
## linedata = alldata

## lineplot = qplot(log(cost), log(wcost),
##   data=linedata,
##   color=as.character(noutput),
##   geom='line',
##   group=as.character(noutput))

## lineplot = lineplot + facet_grid(strat ~ qsize, labeller=lbl.fn)
## lineplot = lineplot +  scale_color_discrete(name='# output cells')
## lineplot = lineplot +  opts(
##   strip.text.y = theme_text(size=6, angle = -90),
##   plot.title = theme_text(size=10))
## lineplot = lineplot +  scale_x_continuous('Provenance Query Cost (sec), log scale')
## lineplot = lineplot +  scale_y_continuous('Operator Overhead (sec), log scale')


## pdf(file='_figs/compare_forward.pdf', width=10, height=10)
## lineplot  %+% linedata[linedata$backward==0,] + opts(title='Query performance vs Runtime Overhead\nForward Queries')
## dev.off()


## pdf(file='_figs/compare_backward.pdf', width=10, height=10)
## lineplot  %+% linedata[linedata$backward==1,] + opts(title='Query performance vs Runtime Overhead\nBackward Queries')
## dev.off()


## paydata = alldata[!is.na(alldata$payload_size),]
## qplot(fanout, disk / 1048576, data=bpaydata, colour=payload_size, geom='line', group=payload_size) + facet_grid(noutput ~ strat, scales='free_y')
## paydata$strat = str_extract(paydata$strat, '^[^_]+_[^_]+')
## qplot(fanout, cost, data=bpaydata[bpaydata$noutput == 100000,], colour=strat, geom='smooth', group=strat) +facet_grid(qsize ~ ., labeller=lbl.fn)












pdf(file='_figs/noutput.pdf', width=10, height=4)
oplot = qplot(noutput, stat, data=sdata, group=strat, colour=strat, shape=strat, linetype=strat, geom='line')
#oplot = oplot + geom_hline(y=1000*1000*4/1048576, alpha=0.7)
oplot = oplot + facet_grid(. ~ fanin, labeller=lbl.fn)
oplot = oplot + scale_x_continuous('% of output cells with provenance')
oplot = oplot + scale_y_continuous('Disk Space Overhead (MB)')
oplot = oplot + opts(
  title='Disk Space Overhead (MB)',
  legend.position='bottom')
oplot = oplot + guides(colour=guide_legend(ncol=4))
oplot
dev.off()


pdf(file='_figs/noutput_overhead.pdf', width=10, height=4)
oplot = qplot(noutput, log(wcost), data=sdata, group=strat, colour=strat, geom='line')
oplot = oplot + facet_grid(. ~ fanin, labeller=lbl.fn)
oplot = oplot + scale_x_continuous('% of output cells with provenance')
oplot = oplot + scale_y_continuous('Runtime Overhead (sec)')
oplot = oplot + opts(
  title='Runtime Overhead',
  legend.position='bottom')
oplot = oplot + guides(colour=guide_legend(ncol=4))
oplot
dev.off()

library('RSQLite')
library('ggplot2')
library(stringr)
require(grid)


theme_set(theme_bw())

lbl.fn = function(var, val) {

  if (var == 'qsize') {
    paste('Query size:', val)
  } else if (var == 'fanout') {
    paste('Fanout:', val)
  } else if (var == 'fanin') {
    paste('Fanin:', val)
  } else {
    print(paste(var, val))
    as.character(val)
  }
}

breaks.fn <- function(x){
    breaks <- c(min(x),median(x),max(x))
    names(breaks) <- attr(breaks,"labels")
    breaks
}

fmt.fn = function(p) {
  p = p + opts(
    legend.text = theme_text(size=15),
    legend.title = theme_text(size=17),
    legend.background = theme_rect(colour = 'white', fill='white', size=0),
    legend.position='bottom',
    #panel.grid.major = theme_line(),
    panel.grid.minor = theme_blank(),
    panel.backgroud = theme_blank(),
    #panel.background = theme_rect(colour = NA),
    strip.background= theme_blank(),
    axis.text.x=theme_text(size=13, colour='#666666'),
    axis.text.y=theme_text(size=13, angle=0, colour='#666666'),
    axis.title.x=theme_text(size=17, lineheight=0.8),
    axis.title.y=theme_text(size=17, angle=90, lineheight=0.8),
    strip.text.x=theme_text(size=15),
    strip.text.y=theme_text(size=15, angle=-90),
    plot.margin = unit(rep(0, 4), "lines"))

  p = p + guides(shape=guide_legend(ncol=4), colour=guide_legend(ncol=4))

  p
}

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#E64341", "#0072B2", "#D55E00", "#CC79A7", "#000000")

drv = dbDriver('SQLite')
con = dbConnect(drv, dbname='./results/micro.jul.9.2012')
alldata = dbGetQuery(con, 'select * from stats, qcosts where stats.id = qcosts.sid')
alldata = alldata[!(alldata$strat %in% c('MANY_KEY_b', 'ONE_MANY_b')),]
#alldata = alldata[!str_detect(alldata$strat, 'PTMAP_ONE'),]
alldata$Strategy = alldata$strat

pidxs = str_detect(alldata$Strategy, 'PTMAP')
alldata$payload_size[alldata$payload_size == 0] = 1
alldata$fanin[pidxs] = alldata$payload_size[pidxs]
alldata$Strategy[str_detect(alldata$Strategy, 'PTMAP_MANY')] =   '<-, PayMany'
alldata$Strategy[str_detect(alldata$Strategy, 'PTMAP_ONE')] =    '<-, PayOne'
alldata$Strategy[str_detect(alldata$Strategy, 'MANY_MANY_b')] =  '<-, FullMany'
alldata$Strategy[str_detect(alldata$Strategy, 'ONE_KEY_b')] =    '<-, FullOne'
alldata$Strategy[str_detect(alldata$Strategy, 'ONE_KEY_f')] =    '->, FullOne'
alldata$Strategy[str_detect(alldata$Strategy, 'Q')] =            'BlackBox' 
levels = c('<-, PayMany',
           '<-, PayOne',
           '<-, FullMany',
           '<-, FullOne',
           '->, FullOne',
           'BlackBox')
alldata$Strategy = factor(alldata$Strategy, levels)


## data = alldata
## data = data[data$fanout == 1 & data$fanin %in% c(1, 100) & data$noutput == 100000,]
## p = qplot(noutput, wcost, data=data, color=Strategy, group=Strategy, geom='line')
## p = p + facet_grid(. ~ fanin, scale='free_y', labeller=lbl.fn)
## p = p + scale_x_continuous('% Output with Lineage')
## p = p + scale_y_continuous('Runtime (sec, log)', limits=c(0, 100))
## p = p + opts(
##   title = 'Runtime',
##   legend.text = theme_text(size=15),
##   legend.title = theme_text(size=15),
##   plot.margin = unit(rep(0, 4), "lines"))
## pdf(file= '_figs/noutput.pdf', width=10, height=3)
## p
## dev.off()




bigdata = alldata
bigdata = bigdata[bigdata$noutput %in% c(100000),]
bigdata = bigdata[!(bigdata$fanout %in% c(2, 5)),]
#bigdata = bigdata[bigdata$Strategy != 'BlackBox',]


data = bigdata[bigdata$fanout %in% c(1, 10), c('fanout', 'Strategy', 'fanin', 'disk', 'wcost')]
dd = data.frame(fanin=c(1, 10))
dd$fanout = 10
dd$Strategy = 'BlackBox'
dd$disk = 0
dd$wcost = 0
data = rbind(data, dd)
d1 = data
d1$stat = d1$disk / 1048576
d1$statname = 'Disk'
d2 = data
d2$stat = d2$wcost
d2$statname = 'Runtime'

data = rbind(d1, d2)


p = qplot(fanin, stat, data=data, color=Strategy, shape=Strategy, group=Strategy, geom='line')
p = p + geom_point()
p = p + scale_shape(solid=F)
p = p + scale_shape_discrete()
p = p + facet_grid(statname ~ fanout, scale='free_y', labeller=lbl.fn)
#p = p + coord_cartesian(ylim=c(0, 40)) 
p = p + scale_color_manual(values=cbbPalette)
p = p + scale_x_continuous('Fanin')
p = p + scale_y_continuous('Runtime (sec)             Disk (MB)\n', breaks=c(0, 10, 20, 30))
p = fmt.fn(p)
pdf(file= '_figs/overhead.pdf', width=10, height=5)
p
dev.off()



## p = qplot(fanin, wcost, data=bigdata, color=Strategy, shape=Strategy, group=Strategy, geom='line')
## p = p + geom_point()
## p = p + scale_shape(solid=F)
## p = p + scale_shape_discrete()
## p = p + facet_grid(. ~ fanout, scale='free_y', labeller=lbl.fn)
## p = p + coord_cartesian(ylim=c(0, 40)) 
## p = p + scale_color_manual(values=cbbPalette)
## p = p + scale_x_continuous('Fanin')
## p = p + scale_y_continuous('Runtime (sec)')#, limits=c(0, 150))
## p = fmt.fn(p)
## pdf(file= '_figs/runtime.pdf', width=10, height=3)
## p
## dev.off()



## p = qplot(fanin, disk / 1048576.0, data=bigdata, color=Strategy, group=Strategy, shape=Strategy, solid=F, geom=c('line', 'point'))
## p = p + geom_hline(y = 1000 * 1000 * 4 / 1048576.0, size=0.3, colour="#555555", linetype=2)
## p = p + facet_grid(. ~ fanout, scale='free_y', labeller=lbl.fn)
## p = p + scale_color_manual(values=cbbPalette)
## p = p + scale_shape_discrete()
## p = p + scale_x_continuous('Fanin')
## p = p + scale_y_continuous('Disk (MB)')
## p = p + coord_cartesian(ylim=c(0, 30)) 
## p = fmt.fn(p)
## pdf(file= '_figs/disk.pdf', width=10, height=3)
## p
## dev.off()




qdata = alldata
qdata = qdata[qdata$fanin %in% c(1, 50, 100) & qdata$qsize == 1000 & qdata$noutput == 100000,]


bdata = qdata[qdata$backward == 1 & str_detect(qdata$Strategy, '<-') & qdata$fanout %in% c(1, 100),]
p = qplot(fanin, cost, data=bdata, group=Strategy, color=Strategy, shape=Strategy, solid=F, geom=c('point','line'))
p = p + facet_grid(.~fanout, labeller=lbl.fn)
p = p + scale_color_manual(values=cbbPalette, name = 'Strategy')
p = p + scale_shape_discrete(name = 'Strategy')
p = p + guides(shape=F)
p = p + scale_x_continuous('Fanin')
p = p + scale_y_continuous('Query Cost (sec)\n', breaks=c(0, .05, .1))#, limits=c(0.0001,1))
p = fmt.fn(p)
pdf(file= '_figs/query_back.pdf', width=10, height=3)
p
dev.off()




fdata = qdata[qdata$backward == 0 & str_detect(qdata$Strategy, '->') & qdata$fanout %in% c(1, 100),]
fdata$fanout = factor(as.character(fdata$fanout), c('1', '50', '100'))
p = qplot(fanin, cost, data=fdata, group=fanout, solid=F, geom=c('point','line'))
p = p + facet_grid(.~fanout, labeller=lbl.fn)
#p = p + scale_color_manual(values=cbbPalette, name = 'Fanout')
#p = p + scale_shape_discrete(name = 'Fanout')
p = p + scale_x_continuous('Fanin')
p = p + scale_y_continuous('Query Cost (sec)\n', breaks=c(0, 0.03, 0.06))
p = fmt.fn(p)
p = p + opts(axis.title.y=theme_text(size=17, angle=90, lineheight=0.8, hjust=0.9))
pdf(file= '_figs/query_forw.pdf', width=10, height=3)
p
dev.off()


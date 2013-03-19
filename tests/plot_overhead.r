library('RSQLite')
library('ggplot2')
library('stringr')
require(grid)
require(viewport)
source('plot_util.r')

theme_set(theme_bw())

args = commandArgs(T)
if (length(args) < 3) {
  print("dude, give me:  runmode, dbpath, outfilepath")
  quit()
}


RUNMODE = as.integer(args[1])
dbName = args[2]
outName = args[3]

drv = dbDriver('SQLite')
con = dbConnect(drv, dbName)



exps = getExperiments(con, RUNMODE)
overhead = getOverhead(con, exps)

if (sum(exps$runtype == 'SubZero10') > 0) {
  levs = c('BlackBox',
    'SubZero1',
    'SubZero10',
    'SubZero20',
    'SubZero50',
    'SubZero100')
  exps$runtype = factor(exps$runtype, levs)
}


merged = merge(exps, overhead, on='id')

# for NSF plots
## merged = merged[merged$runtype %in% c('BlackBox', 'FullOne', 'SubZero'),]
## merged$runtype[merged$runtype == 'FullOne'] = 'Cell Level'
## merged

d1 = merged
d1$stat = d1$disk
d1$statname = 'Disk Cost'
d2 = merged
d2$stat = d2$opcost
d2$statname = 'Runtime'
compdata = rbind(d1, d2)
compdata$Strategy = compdata$runtype

p = ggplot(compdata, aes(x=Strategy, y=stat, fill=Strategy, label=as.integer(compdata$stat)))
p = p + geom_bar(position='dodge', stat='identity', aes(width=0.5))
p = p + geom_text(aes(y=max(stat)*0.9))
p = p + facet_grid(statname ~ .)
#p = p + scale_y_log10('Runtime (sec)     Disk (MB)')#, limits=c(1, 3)) 
p = p + scale_y_continuous('Runtime       Disk\n  (sec)        (MB)\n')
#p = p + coord_cartesian(ylim=c(0, 77)) 
#p = p + scale_fill_brewer(palette='OrRd')
p = p + scale_fill_manual(values=cbbPalette)
p = p + scale_x_discrete('Storage Strategies')
p = p + opts(plot.margin = unit(rep(0.1, 4), "lines"))
p = fmt.fn(p)
p = p + opts(axis.title.y=theme_text(size=17, angle=90, lineheight=0.8, hjust=0.5))
pdf(file=outName, width=10, height=3.5)
p
dev.off()

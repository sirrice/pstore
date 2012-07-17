library('RSQLite')
library('ggplot2')
library('stringr')
require(grid)
require(viewport)
source('plot_util.r')

theme_set(theme_bw())

args = commandArgs(T)
if (length(args) < 2) {
  print("dude, give me:  runmode, outfilepath")
  quit()
}


RUNMODE = as.integer(args[1])
outName = args[2]


drv = dbDriver('SQLite')
con = dbConnect(drv, dbname='./results/lsstfull.db.nov.24.2011.fastforward')


exps = getExperiments(con, RUNMODE)

r = getQueryCosts(con, exps)
compdata = merge(exps, r)

fcon = dbConnect(drv, dbname='./results/lsstfull.db.nov.24.2011.slowforward')
fexps = getExperiments(fcon, RUNMODE)
fr = getQueryCosts(fcon, fexps)
fr = fr[fr$forward == 1,]
fr$name = paste(fr$name, "Slow")
fexps = fexps[fexps$id %in% fr$id,]
fcompdata = merge(fexps, fr, on='id')

data = rbind(compdata, fcompdata)
data$Strategy = data$runtype
data$Query = data$name

p = qplot(Strategy, cost, data=data, fill=Query, geom='bar', position='dodge')
p = p + scale_y_log10('Query Cost\n(sec, log)\n')
p = p + scale_fill_manual(values=cbbPalette)
p = p + scale_x_discrete('Storage Strategies')
p = fmt.fn(p)
#p = p + opts(axis.title.y=theme_text(size=17, angle=90, lineheight=0.8, hjust=0.9))

pdf(file=outName, width=10, height=3)
p
dev.off()




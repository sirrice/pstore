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
con = dbConnect(drv, dbname=dbName)


exps = getExperiments(con, RUNMODE)

r = getQueryCosts(con, exps)
data = merge(exps, r)


if (sum(exps$runtype == 'SubZero10') > 0) {
  levs = c('BlackBox',
    'SubZero1',
    'SubZero10',
    'SubZero20',
    'SubZero50',
    'SubZero100')
  exps$runtype = factor(exps$runtype, levs)
}
colnames(data)
data$Strategy = data$runtype
data$Query = data$name

p = ggplot(data, aes(x=Strategy, y=cost, fill=Query))
p = p + geom_bar(position='dodge', stat='identity', aes(width=0.7))
p = p + scale_y_log10('Query Cost\n(sec, log)\n')
p = p + scale_fill_manual(values=cbbPalette)
p = p + scale_x_discrete('Storage Strategies')
p = fmt.fn(p)
#p = p + opts(axis.title.y=theme_text(size=17, angle=90, lineheight=0.8, hjust=0.9))
pdf(file=outName, width=10, height=3)
p
dev.off()




library('RSQLite')
library('ggplot2')


drv = dbDriver('SQLite')
con = dbConnect(drv, dbname='./results/micro.jul.9.2012')

dbListTables(con)

alldata = dbGetQuery(con, 'select * from stats, qcosts where stats.id = qcosts.sid')

colnames(alldata)



lineplot = qplot(log(cost), log(wcost),
  data=alldata,
  color=as.character(noutput),
  facets=strat ~ qsize,
  geom='line',
  group=as.character(noutput))

lineplot = lineplot +  scale_color_discrete(name='# output cells')

lineplot = lineplot +  opts(
  strip.text.y = theme_text(size=6, angle = -90),
  plot.title = theme_text(size=10))

lineplot = lineplot +  scale_x_continuous('Provenance Query Cost (sec), log scale')
lineplot = lineplot +  scale_y_continuous('Operator Overhead (sec), log scale')


pdf(file='forward_plot.pdf', width=10, height=10)
lineplot  %+% alldata[alldata$backward==0,] + opts(title='Query performance vs Runtime Overhead\nForward Queries')
dev.off()


pdf(file='backward_plot.pdf', width=10, height=10)
lineplot  %+% alldata[alldata$backward==1,] + opts(title='Query performance vs Runtime Overhead\nBackward Queries')
dev.off()




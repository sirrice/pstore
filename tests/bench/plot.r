library(RSQLite)
library(ggplot2)
library(sqldf)

fname = commandArgs(trailingOnly=T)[1]
db = dbConnect(SQLite(), dbname = fname)
stats = sqldf(conn=db, "select * from stats")
qcosts = sqldf(conn=db, "select * from qcosts")

stats$disk = stats$disk / 1048576

data = stats

p = ggplot(data, aes(x=fanin, y=disk, group=strat, color=strat)) + facet_wrap(noutput~fanout) + geom_line() + coord_cartesian(ylim=c(0, 30))
ggsave(filename=paste(fname, ".plot.pdf", sep=''), plot=p)
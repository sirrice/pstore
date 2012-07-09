library('RSQLite')
library('ggplot2')


drv = dbDriver('SQLite')
#con = dbConnect(drv, dbname='./results/pstore_microbench.db.jul.8.201')
con = dbConnect(drv, dbname='/tmp/db.db')

dbListTables(con)

alldata = dbGetQuery(con, 'select * from stats, qcosts where stats.id = qcosts.sid')

colnames(alldata)



qplot(cost, wcost, data=alldata, color=strat, shape=as.character(backward))

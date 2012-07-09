library('RSQLite')
# library('ggplot')


drv = dbDriver('SQLite')
con = dbConnect(drv, dbname='./results/pstore_microbench.db.jul.8.201')

dbListTables(con)

alldata = dbGetQuery(con, 'select * from stats, qcosts where stats.id = qcosts.sid')

colnames(alldata)

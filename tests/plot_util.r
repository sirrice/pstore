cbbPalette <- c("#E69F00", "#56B4E9", "#00af84", "#E64341", "#0072B2", "#D55E00", "#CC79A7", "#000000")

getExperiments = function(con, runmode) {
  q = sprintf("select rowid, runmode, runtype, width, height, diskconstraint, runconstraint
                from exec where runmode = %d and
                runtype not in ('noop', 'noop_model', 'stats', 'stats_model', 'opt', 'noop_m', 'stats_m')
                order by diskconstraint", runmode)
  res = dbGetQuery(con, q)
  res$runtype = toupper(res$runtype)
  
  replaces = data.frame(
    k=c('KEY', 'PTR_', '_B', '_F', 'PTMAP_', 'OPT_', '.0_[0-9.]+$', 'MANUAL$'),
    v=c('REF', '3_', '_b', '_f', '2_', 'SubZero', '', ''))
  for (rid in rownames(replaces)) {
    res$runtype = sub(replaces[rid, 'k'], replaces[rid, 'v'], res$runtype)
  }

  replaces = data.frame(
    k = c('Q_OPT', 'F_B', 'OPT_MANUAL', 'QUERY', 'ONE_REF', 'MANY_REF', 'MANY_MANY', 'QUERYOPT', '3_ONE_REF_B', '3_MANY_MANY_B', '3_ONE_REF_F', '3_F_B', '2_ONE_REF_B', '2_MANY_MANY_B', '2_F_B','MANY_MANY', 'OPT_', 'Q_OPT'),
    v=c('BlackBox', 'ONE_REF_b,f', 'SubZero', 'BlackBox', 'FullOne', 'FullMany', 'FullMany', 'BlackBoxOpt', 'FullOne', 'FullMany', 'FullForw', 'FullBoth', 'PayOne', 'PayMany', 'PayBoth', 'FullMany', 'SubZero', 'BlackBox'))
  
  for (rid in rownames(replaces)) {
    #res$runtype = sub(replaces[rid, 'k'], replaces[rid, 'v'], toupper(res$runtype))
    idxs = tolower(res$runtype) == tolower(replaces[rid, 'k'])
    if (sum(idxs) > 0) {
      print(paste(res$runtype[idxs][1], '-->', replaces[rid, 'v'], ' on ', replaces[rid, 'k']))
      res$runtype[idxs] = as.character(replaces[rid, 'v'])
    }
  }

  res
}

vdbGetQuery = function(con, q, df, cols) {
  ret = data.frame()
  for (rowname in rownames(df)) {
    ret = rbind(ret, dbGetQuery(con, sprintf(q, unlist(df[rowname,cols]))))
  }
  ret  
}


getOverhead = function(con, exps) {
  q = "SELECT exec.id, sum(opcost) as opcost, sum(disk) as disk, 0 as idx
    FROM pstore_overhead as po, workflow_run as wr, exec
    WHERE po.wid = wr.rowid and wr.eid = exec.rowid and exec.rowid = %d and strat != 's:Func' "
  overheads = vdbGetQuery(con, q, exps, 'id')

  
  q = "SELECT sum(wr.opcost) as opcost, 0 as disk, 0 as idx
    FROM workflow_run as wr, exec as e1, exec as e2
    WHERE e1.rowid = wr.eid and e1.runtype = 'noop' and
    strat != 's:Func' and e2.rowid = %d and e1.runmode = e2.runmode;"
  baselines = vdbGetQuery(con, q, exps, 'id')

  baselines$disk = exps$width * exps$height * 8.00002 * 2 / 1048576
  overheads$disk = overheads$disk / 1048576. + baselines$disk#8.6
  overheads  
}




getQueryCosts = function (con, exps) {

  q = "SELECT pq.eid as id, pq.rowid, pq.forward, pp.opid, pq.cost as cost
 FROM pq, pq_path as pp
 WHERE pp.pqid = pq.rowid and pq.eid = %d
 ORDER BY pq.forward, pq.rowid, pp.idx;"

  r = vdbGetQuery(con, q, exps, 'id')
            
  r = aggregate(r, list(r$rowid), rbind)
  r$rowid = sapply(r$rowid, function(x)x[1])
  r$forward = sapply(r$forward, function(x)x[1])
  r$id = sapply(r$id, function(x)x[1])
  r$opid = sapply(r$opid, function(l)paste(unlist(l), collapse=','))
  r$cost = sapply(r$cost, mean)
  r = subset(r, select=-c(Group.1))

  # construct names 
  allpaths = aggregate(r, list(r$opid), rbind)
  allpaths$forward = sapply(allpaths$forward, function(x)x[1])
  allpaths$opid = sapply(allpaths$opid, function(x)x[1])

  # forward path names
  fidxs = allpaths$forward == 1
  allpaths$name[fidxs] = paste('FQ ', 1:length(fidxs) - 1, sep='')
  
  # backward path names
  bidxs = allpaths$forward == 0
  allpaths$name[bidxs] = paste('BQ ', 1:length(bidxs) - 1, sep='')
  
  #allpaths$name = paste(ifelse(allpaths$forward == 1, 'F', 'B'), 'Q ', rownames(allpaths), sep='')
  rownames(allpaths) = allpaths$opid

  r$name = allpaths[r$opid, 'name']
  r
}




fmt.fn = function(p) {
  p = p + opts(
    legend.text = theme_text(size=15),
    legend.title = theme_text(size=17),
    legend.background = theme_rect(colour = 'white', fill='white', size=0),
    legend.position='bottom',
    legend.key.width = unit(0.7, "cm"),
    #panel.grid.major = theme_line(),
    panel.grid.minor = theme_blank(),
    panel.backgroud = theme_blank(),
    #panel.background = theme_rect(colour = NA),
    strip.background= theme_blank(),
    axis.text.x=theme_text(size=14, colour='#666666'),
    axis.text.y=theme_text(size=14, angle=0, colour='#666666'),
    axis.title.x=theme_text(size=17, lineheight=0.8),
    axis.title.y=theme_text(size=17, angle=90, lineheight=0.8),
    strip.text.x=theme_text(size=15),
    strip.text.y=theme_text(size=15, angle=-90)
    #plot.margin = unit(rep(0, 4), "lines")
    )

  p = p + guides(shape=guide_legend(ncol=4), colour=guide_legend(ncol=4), fill=guide_legend(ncol=4))

  p
}

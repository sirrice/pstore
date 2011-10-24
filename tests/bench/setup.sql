
drop view if exists ptr_plot_abs ;
create view ptr_plot_abs as 
select stats.rid, strat, fanin, fanout, density, nqs, noutput,
       disk, overhead, sercost, writecost, fcost, bcost, round(avgoc.min::numeric,4) as opcost
from runs, stats, (select s2.rid, min(s2.opcost) from stats as s2 group by rid) as avgoc 
where avgoc.rid = stats.rid and runs.id = stats.rid and runs.notes = 'ptr'
order by stats.rid;



drop view if exists model_plot_abs ;
create view model_plot_abs as 
select stats.rid, strat, fanin, fanout, density, nqs, noutput,
       disk, overhead, sercost, writecost, fcost, bcost, round(avgoc.min::numeric,4) as opcost
from runs, stats, (select s2.rid, min(s2.opcost) from stats as s2 group by rid) as avgoc 
where avgoc.rid = stats.rid and runs.id = stats.rid and runs.notes = 'ptr_model'
order by stats.rid;


-- drop view if exists ptr_plot_perc ;
-- create view ptr_plot_perc as 
-- select stats.rid, strat, fanin, fanout, density, nqs, noutput,
--        disk/arrsize::float as disk, round((overhead/avgoc.min)::numeric,4) as overhead, 
--        round((sercost/opcost)::numeric, 4) as sercost,
--        round((writecost/opcost)::numeric, 4) as writecost,
--        fcost, bcost, round(avgoc.avg::numeric,4) as opcost, arrsize
-- from runs, stats, (select s2.rid, min(s2.opcost) from stats as s2 group by rid) as avgoc 
-- where avgoc.rid = stats.rid and runs.id = stats.rid and runs.notes = 'ptr'
-- order by stats.rid;


drop view if exists box_plot_abs ;
create view box_plot_abs as 
select stats.rid, strat, fanin, fanout, density, nqs, noutput,
       disk, overhead, sercost, writecost, fcost, bcost, round(opcost::numeric,3) as opcost
from runs, stats
where runs.id = stats.rid and runs.notes = 'box'
order by stats.rid;

-- drop view if exists box_plot_perc ;
-- create view box_plot_perc as 
-- select stats.rid, strat, fanin, fanout, density, nqs, noutput,
--        disk/arrsize::float as disk, round((overhead/opcost)::numeric,4) as overhead, 
--        round((sercost/opcost)::numeric, 4) as sercost,
--        round((writecost/opcost)::numeric, 4) as writecost,
--        round((fcost/opcost)::numeric,4) as fcost, round((bcost/opcost)::numeric,4) as bcost, 
--        round((opcost::numeric),3) as opcost
-- from runs, stats
-- where runs.id = stats.rid and runs.notes = 'box'
-- order by stats.rid;



drop view if exists box_model_plot_abs ;
create view box_model_plot_abs as 
select stats.rid, strat, fanin, fanout, density, nqs, noutput,
       disk, overhead, sercost, writecost, fcost, bcost, round(opcost::numeric,3) as opcost
from runs, stats
where runs.id = stats.rid and runs.notes = 'box_model'
order by stats.rid;




-- slope of serialization cost vs fanin
drop table if exists sercostvsfanin ;
create table sercostvsfanin as
select rid, strat, fanout, density, nqs, noutput, opcost, 
       max(sercost) as maxser, min(sercost) as minser, 
       max(fanin) as maxfanin, min(fanin) as minfanin, 
       (max(sercost)-min(sercost)) / (max(fanin)-min(fanin)+1) as serslope 
from ptr_plot_abs 
group by rid, strat, fanout, density, nqs, noutput, opcost 
having max(fanin)-min(fanin) > 0 and max(sercost) > 0 and min(sercost) > 0 
order by strat, fanout;

select rid, strat, density, nqs, noutput, opcost, 
       stddev(serslope) as std, avg(serslope) as avg, min(serslope) as min,
       (max(serslope)-min(serslope))/(max(fanout)-min(fanout)) as deltaslope
from sercostvsfanin 
where rid = 16
group by rid, strat, density, nqs, noutput, opcost;



-- slope of write cost vs fanin
drop table if exists writecostvsfanin ;
create table writecostvsfanin as
select rid, strat, fanout, density, nqs, noutput, opcost, 
       max(writecost) as maxwrite, min(writecost) as minwrite, 
       max(fanin) as maxfanin, min(fanin) as minfanin, 
       (max(writecost)-min(writecost)) / (max(fanin)-min(fanin)+1) as writeslope 
from ptr_plot_abs 
group by rid, strat, fanout, density, nqs, noutput, opcost 
having max(fanin)-min(fanin) > 0 and max(writecost) > 0 and min(writecost) > 0 
order by strat, fanout;

select rid, strat, density, nqs, noutput, opcost, 
       stddev(writeslope) as std, avg(writeslope) as avg, min(writeslope) as min,
       (max(writeslope)-min(writeslope))/(max(fanout)-min(fanout)) as deltaslope
from writecostvsfanin 
where rid = 16
group by rid, strat, density, nqs, noutput, opcost;


-- slope of overhead vs fanin
drop table if exists overheadvsfanin ;
create table overheadvsfanin as
select rid, strat, fanout, density, nqs, noutput, opcost, 
       max(overhead) as maxovh, min(overhead) as minovh, 
       max(fanin) as maxfanin, min(fanin) as minfanin, 
       (max(overhead)-min(overhead)) / (max(fanin)-min(fanin)+1) as ovhslope 
from ptr_plot_abs 
group by rid, strat, fanout, density, nqs, noutput, opcost 
having max(fanin)-min(fanin) > 0 and max(overhead) > 0 and min(overhead) > 0 
order by strat, fanout;

select rid, strat, density, nqs, noutput, opcost, 
       stddev(ovhslope) as std, avg(ovhslope) as avg, min(ovhslope) as min,
       (max(ovhslope)-min(ovhslope))/(max(fanout)-min(fanout)) as deltaslope
from overheadvsfanin 
where rid = 16
group by rid, strat, density, nqs, noutput, opcost;



-- slope of disk cost vs fanin
drop table if exists diskvsfanin ;
create table diskvsfanin as
select rid, strat, fanout, density, nqs, noutput, opcost, 
       max(disk) as maxdisk, min(disk) as mindisk, 
       max(fanin) as maxfanin, min(fanin) as minfanin, 
       (max(disk)-min(disk)) / (max(fanin)-min(fanin)+1) as diskslope 
from ptr_plot_abs 
group by rid, strat, fanout, density, nqs, noutput, opcost 
having max(fanin)-min(fanin) > 0 and max(disk) > 0 and min(disk) > 0 
order by strat, fanout;

select rid, strat, density, nqs, noutput, opcost, 
       stddev(diskslope) as std, avg(diskslope) as avg, min(diskslope) as min,
       (max(diskslope)-min(diskslope))/(max(fanout)-min(fanout)) as deltaslope
from diskvsfanin 
where rid = 16
group by rid, strat, density, nqs, noutput, opcost;
#!/bin/bash

r -f plot_overhead.r --args 1 results/lsstfull.db.nov.24.2011.fastforward figs/lsst/overhead.pdf
r -f plot_lsst.r --args 1 figs/lsst/cost.pdf


r -f plot_overhead.r --args 100 results/pablostats_load.db.nov.21.2011 figs/genomics/overhead_noopt.pdf 
r -f plot_pablo.r --args 100 results/pablostats_static.db.nov.20.2011 figs/genomics/cost_static.pdf
r -f plot_pablo.r --args 100 results/pablostats_dynamic.db.nov.20.2011 figs/genomics/cost_dynamic.pdf

r -f plot_overhead.r --args 100 results/pablostats_100_opt.db.nov.30.2011 figs/genomics/overhead_opt.pdf
r -f plot_pablo.r --args 100 results/pablostats_100_opt.db.nov.30.2011 figs/genomics/cost_opt.pdf

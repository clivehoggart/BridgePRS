library(BEDMatrix)
library(MASS)
library(parallel)
library(data.table)
library("optparse")
options(stringsAsFactors=FALSE)

#option_list = list(
#    make_option(c("--fpath"), type="character", default=NULL,help="Function File Path", metavar="character"),
#    make_option(c("--by.chr.sumstats"), type="character", default=0,help="Logical indicating if sumstats files are split by chr",metavar="character"))
#
#opt_parser = OptionParser(option_list=option_list);
#opt = parse_args(opt_parser)



#print(opt)

args <- commandArgs(trailingOnly=TRUE) 

logfile <- args[2] 

packages <- c('ggplot2','dplyr','lattice') 

#packages <- c('dplyr','ggplot','BEDMatrix','MASS','parallel','data.table','optparse') 
packages <- c('BEDMatrix','MASS','parallel','data.table','R.utils','optparse','boot','doMC','glmnet') 

needed_packages = setdiff(packages, rownames(installed.packages())) 
#needed_packages <- c('boot','MASS') 
#needed_packages <- c('boot') 

if (args[1] == 'CHECK') 
    { 
        write(needed_packages, logfile, sep='\n', append=TRUE) 
} else {
    install.packages(needed_packages, repos = "http://cran.us.r-project.org")
}
quit() 

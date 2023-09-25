args <- commandArgs(trailingOnly=TRUE) 




logfile <- args[2] 

packages <- c('ggplot2','dplyr','lattice') 

#packages <- c('dplyr','ggplot','BEDMatrix','MASS','parallel','data.table','optparse') 
packages <- c('BEDMatrix','MASS','parallel','data.table','R.utils','optparse','boot','doMC','glmnet') 

packages <- c('BEDMatrix','MASS','parallel','data.table','R.utils','optparse','boot','doMC','glmnet') 


needed_packages = setdiff(packages, rownames(installed.packages())) 


if (args[1] == 'CHECK') 
    { 
        write(needed_packages, logfile, sep='\n', append=TRUE) 
} else {
    install.packages(needed_packages, repos = "http://cran.us.r-project.org")
}
quit() 

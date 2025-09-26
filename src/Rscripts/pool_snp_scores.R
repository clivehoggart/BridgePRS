library(data.table)
library("optparse")

option_list = list(
    make_option(c("--n.folds"), type="numeric", default=0,
                help="cv folds", metavar="numeric"),
    make_option("--workdir", type="character", default=NULL,
                help="file path", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print(opt)

model <- c('best','weighted')

for( i in 1:2 ){
    snps <- vector()
    tmp <- list()
    for( iter in 1:opt$n.folds ){
        tmp[[iter]] <- fread(paste0(opt$workdir,'/fold',iter,'/snp_weights_',model[i],'_model.dat'))
        snps <- union( snps, tmp[[iter]]$SNP.id )
    }
    ensembl.model <- as.data.frame(matrix( ncol=4, nrow=length(snps), data=0 ))
    ensembl.model[,1] <- snps
    for( iter in 1:opt$n.folds ){
        ptr <- match( tmp[[iter]]$SNP.id, snps )
        ensembl.model[ptr,2:3] <- c( tmp[[iter]]$ALLELE1, tmp[[iter]]$ALLELE0 )
        ensembl.model[ptr,4] <- ensembl.model[ptr,4] + tmp[[iter]][,4]
    }
    ptr.use <- which(ensembl.model[,4]!=0)
    ensembl.model <- ensembl.model[ptr.use,]
    colnames(ensembl.model) <- c('SNP.id', 'ALLELE1', 'ALLELE0', 'beta' )
    fwrite( ensembl.model,
           paste0(opt$workdir,'/snp_weights_',model[i],'_model.dat'),
           row.names=FALSE, col.names=TRUE, sep=" " )
}

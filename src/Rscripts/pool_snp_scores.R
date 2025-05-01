library(data.table)

option_list = list(
    make_option("--workdir", type="character", default=NULL,
                help="file path", metavar="character")
)

dir <- 'out/preds/'
pop <- 'EAS'
model <- c('best','weighted')

for( i in 1:2 ){
    snps <- vector()
    tmp <- list()
    for( iter in 1:5 ){
        tmp[[iter]] <- fread(paste0(opt$workdir,'_snp_weights_',model[i],'_model_',iter,'.dat'))
        snps <- union( snps, tmp[[iter]]$V1 )
    }
    ensembl.model <- as.data.frame(matrix( ncol=4, nrow=length(snps), data=0 ))
    ensembl.model[,1] <- snps
    for( iter in 1:5 ){
        ptr <- match( tmp[[iter]]$V1, snps )
        ensembl.model[ptr,2:3] <- c( tmp[[iter]]$V2, tmp[[iter]]$V3 )
        ensembl.model[ptr,4] <- ensembl.model[ptr,4] + tmp[[iter]]$V4
    }
    ptr.use <- which(ensembl.model[,4]!=0)
    ensembl.model <- ensembl.model[ptr.use,]
    write.table( ensembl.model,
                paste0(opt$workdir,'_snp_weights_',model[i],'_model_0.dat'),
                row.names=FALSE, col.names=FALSE, quote=FALSE )
}

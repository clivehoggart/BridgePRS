library(BEDMatrix)
library(MASS)
library(parallel)
library(data.table)
library("optparse")
options(stringsAsFactors=FALSE)

# toupper

option_list = list(
    make_option(c("--fpath"), type="character", default=NULL,help="Function File Path", metavar="character"),
    make_option(c("--workdir"), type="character",
                help="Top level working dir", metavar="character"),
    make_option(c("--bfile"), type="character",
                help="Plink file to estimate LD", metavar="character"),
    make_option(c("--by.chr"), type="numeric", default=1,
                help="Logical indicating if bed files are split by chr",
                metavar="character"),
    make_option(c("--ld.ids"), type="character",
                help="File listing IDs to use for estimating LD", metavar="character"),
    make_option(c("--sumstats"), type="character",
                help="Sum stats file", metavar="character"),
    make_option(c("--by.chr.sumstats"), type="character", default=0,
                help="Logical indicating if sumstats files are split by chr",
                metavar="character"),
    make_option(c("--sumstats.snpID"), type="character", default="SNP",
                help="SNP column name", metavar="character"),
    make_option(c("--sumstats.betaID"), type="character", default="BETA",
                help="Beta column name", metavar="character"),
    make_option(c("--sumstats.allele0ID"), type="character", default="ALLELE0",
                help="Allele 0 column name", metavar="character"),
    make_option(c("--sumstats.allele1ID"), type="character", default="ALLELE1",
                help="Allele 1 column name", metavar="character"),
    make_option(c("--sumstats.P"), type="character", default="ALLELE1",
                help="P-value column name", metavar="character"),
    make_option(c("--N.pop"), type="numeric", default=0,
                help="GWAS sample size", metavar="numeric"),
    make_option(c("--prop.train"), type="numeric", default=0,
                help="GWAS training sample proportion", metavar="numeric"),
    make_option(c("--prop.test"), type="numeric", default=0,
                help="Testing sample proportion", metavar="numeric"),
    make_option(c("--strand.check"), type="numeric", default=0,
                help="Keep only non-ambiguous SNPs", metavar="numeric"),
    make_option(c("--n.cores"), type="numeric", default=1,
                help="Number of processors for mclapply to use", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print(opt)

source(opt$fpath)

logfile <- paste0(opt$workdir,"sumstat_subset.log")
tmp <- t(data.frame(opt))
rownames(tmp) <- names(opt)
write.table(tmp,file=logfile,quote=FALSE,col.names=FALSE)

if( opt$by.chr.sumstats==0 ){
    sumstats <- fread( opt$sumstats, data.table=FALSE )
    snp.ptr <- which( colnames(sumstats)==opt$sumstats.snpID )
    allele1.ptr <- which( colnames(sumstats)==opt$sumstats.allele1ID )
    allele0.ptr <- which( colnames(sumstats)==opt$sumstats.allele0ID )
    beta.ptr <- which( colnames(sumstats)==opt$sumstats.betaID )

    sumstats <- sumstats[,c( snp.ptr, allele1.ptr, allele0.ptr, beta.ptr)]
    colnames(sumstats) <- c('SNP','ALLELE1','ALLELE0','BETA')
    sumstats$ALLELE1 <- toupper(sumstats$ALLELE1)
    sumstats$ALLELE0 <- toupper(sumstats$ALLELE0)
    if( opt$strand.check ){
        ptr.use <- which( (sumstats$ALLELE1=="A" & sumstats$ALLELE0=="C") |
                          (sumstats$ALLELE1=="A" & sumstats$ALLELE0=="G") |
                          (sumstats$ALLELE1=="C" & sumstats$ALLELE0=="A") |
                          (sumstats$ALLELE1=="C" & sumstats$ALLELE0=="T") |
                          (sumstats$ALLELE1=="G" & sumstats$ALLELE0=="A") |
                          (sumstats$ALLELE1=="G" & sumstats$ALLELE0=="T") |
                          (sumstats$ALLELE1=="T" & sumstats$ALLELE0=="C") |
                          (sumstats$ALLELE1=="T" & sumstats$ALLELE0=="G") )
        sumstats <- sumstats[ptr.use,]
    }
    sumstats$BETA <- as.numeric(sumstats$BETA)
    sumstats <- sumstats[ !is.na(sumstats$BETA), ]
}

ld.ids <- as.character(read.table(opt$ld.ids)[,2])

if( opt$by.chr==0 ){
    ptr.bed <- BEDMatrix( opt$bfile, simple_names=TRUE )
    bim <- fread( paste(opt$bfile,'.bim',sep='' ) )
    ld.ids <- intersect( ld.ids, attributes(ptr.bed)[[3]][[1]] )
}

for( chr in 1:22 ){
    if( opt$by.chr.sumstats!=0 ){
        sumfile  <- paste(opt$sumstats,chr,opt$by.chr.sumstats,sep='')
        if(!file.exists(sumfile)) {
            next
            }
        sumstats <- fread( paste(opt$sumstats,chr,opt$by.chr.sumstats,sep=''), data.table=FALSE )

        snp.ptr <- which( colnames(sumstats)==opt$sumstats.snpID )
        allele1.ptr <- which( colnames(sumstats)==opt$sumstats.allele1ID )
        allele0.ptr <- which( colnames(sumstats)==opt$sumstats.allele0ID )
        beta.ptr <- which( colnames(sumstats)==opt$sumstats.betaID )
        p.ptr <- which( colnames(sumstats)==opt$sumstats.P )

        sumstats <- sumstats[,c( snp.ptr, allele1.ptr, allele0.ptr, beta.ptr, p.ptr )]

        colnames(sumstats) <- c('SNP','ALLELE1','ALLELE0','BETA','P')
        sumstats$ALLELE1 <- toupper(sumstats$ALLELE1)
        sumstats$ALLELE0 <- toupper(sumstats$ALLELE0)
        if( opt$strand.check ){
            ptr.use <- which( sumstats$ALLELE1=="A" & sumstats$ALLELE0=="C" |
                              sumstats$ALLELE1=="A" & sumstats$ALLELE0=="G" |
                              sumstats$ALLELE1=="C" & sumstats$ALLELE0=="A" |
                              sumstats$ALLELE1=="C" & sumstats$ALLELE0=="T" |
                              sumstats$ALLELE1=="G" & sumstats$ALLELE0=="A" |
                              sumstats$ALLELE1=="G" & sumstats$ALLELE0=="T" |
                              sumstats$ALLELE1=="T" & sumstats$ALLELE0=="C" |
                              sumstats$ALLELE1=="T" & sumstats$ALLELE0=="G" )
            sumstats <- sumstats[ptr.use,]
        }
        sumstats$BETA <- as.numeric(sumstats$BETA)
        sumstats <- sumstats[ !is.na(sumstats$BETA), ]
    }
    if( opt$by.chr==1 ){
        ptr.bed <- BEDMatrix( paste(opt$bfile,chr,sep=''), simple_names=TRUE )
        bim <- fread( paste(opt$bfile,chr,'.bim',sep='' ) )
        ld.ids <- intersect( ld.ids, attributes(ptr.bed)[[3]][[1]] )
    }
    infile <- paste(opt$workdir,'/blocks/chr',chr,'.blocks.det.gz',sep='')
    blocks <- fread(infile,header=TRUE,stringsAsFactors=FALSE)

    sumstats.b <- mclapply( 1:nrow(blocks),
                     function(i){
                         sumstat.subset( block.i=blocks[i,],
                                        sumstats=sumstats, ld.ids=ld.ids,
                                        X.bed=ptr.bed, bim=bim,
                                        n.all=opt$N.pop,
                                        n.prop=c(opt$prop.train, opt$prop.test),
                                        strand.check=opt$strand.check )},
                     mc.cores=as.numeric(opt$n.cores) )
    sumstats.1 <- as.data.frame(matrix(nrow=nrow(sumstats),ncol=7))
    colnames(sumstats.1) <- c('SNP','ALLELE1','ALLELE0','BETA','SE','P','XtY.2','XtY.3')
    for( i in 1:nrow(blocks) ){
        ptr <- match( sumstats.b[[i]]$SNP, sumstats$SNP )
        sumstats.1[ptr,] <- sumstats.b[[i]]
    }
    ptr.fill <- which(is.na(sumstats.1$BETA))
    sumstats.b <- mclapply( ptr.fill,
                     function(i){
                         sumstat.subset( snp=sumstats$SNP[i],
                                        sumstats=sumstats, ld.ids=ld.ids,
                                        X.bed=ptr.bed, bim=bim,
                                        n.all=opt$N.pop,
                                        n.prop=c( opt$prop.train, opt$prop.test ),
                                        strand.check=opt$strand.check )},
                     mc.cores=as.numeric(opt$n.cores) )
    for( i in 1:length(ptr.fill) ){
        ptr <- match( sumstats.b[[i]]$SNP, sumstats$SNP )
        sumstats.1[ptr,] <- sumstats.b[[i]]
    }
    sumstats.1 <- sumstats.1[!is.na(sumstats.1$BETA),]
    outfile <- paste(opt$workdir,'/sumstat_subset/chr',chr,'.dat.gz',sep='')
    fwrite( sumstats.1, outfile, sep=" " )
}
if( !is.null(warnings()) ){
    print(warnings())
}

#for( i in ptr.fill ){
#    tmp <- sumstat.subset( snp=sumstats$SNP[i],
#                          sumstats=sumstats, ld.ids=ld.ids,
#                          X.bed=ptr.bed, bim=bim,
#                          n.all=opt$N.pop, n.gwas=0.7*opt$N.pop,
#                          strand.check=opt$strand.check )
#}
#i=i
#tmp <- sumstat.subset( block.i=blocks[i,],
#                      sumstats=sumstats, ld.ids=ld.ids,
#                      X.bed=ptr.bed, bim=bim,
#                      n.all=opt$N.pop, n.gwas=0.7*opt$N.pop,
#                      strand.check=opt$strand.check )

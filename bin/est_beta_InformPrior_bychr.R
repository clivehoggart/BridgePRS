library(BEDMatrix)
library(MASS)
library(parallel)
library(data.table)
library("optparse")
source('~/BridgePRS/bin/functions.R')
options(stringsAsFactors=FALSE)

option_list = list(
    make_option(c("-c", "--clump.stem"), type="character",
                help="Clump stem", metavar="character"),
    make_option(c("--prior"), type="character",
                help="Prior stem", metavar="character"),
    make_option(c("-s", "--sumstats"), type="character",
                help="Sum stats file", metavar="character"),
    make_option(c("--ld.ids"), type="character",
                help="File listing IDs to use for estimating LD", metavar="character"),
    make_option(c("-b", "--beta.stem"), type="character",
                help="Out file stem", metavar="character"),
    make_option(c("--w.prior"), type="character",
                help="Comma delimited prior weights", metavar="character"),
    make_option(c("--precision"), type="numeric", default=0,
                help="Logicial to calculate precision matrix of each clump",
               metavar="character"),
    make_option(c("--n.cores"), type="numeric", default=1,
                help="Number of processors for mclapply to use", metavar="character"),
    make_option(c("--bfile"), type="character",
                help="Plink file with test & train2 genotype data", metavar="character"),
    make_option(c("--sumstats.snpID"), type="character", default="SNP",
                help="SNP column name", metavar="character"),
    make_option(c("--sumstats.betaID"), type="character", default="BETA",
                help="Beta column name", metavar="character"),
    make_option(c("--sumstats.allele0ID"), type="character", default="ALLELE0",
                help="Allele 0 column name", metavar="character"),
    make_option(c("--sumstats.allele1ID"), type="character", default="ALLELE1",
                help="Allele 1 column name", metavar="character"),
    make_option(c("--sumstats.nID"), type="character", default="OBS",
                help="N. obs column name", metavar="character"),
    make_option(c("--sumstats.seID"), type="character", default="SE",
                help="SE column name", metavar="character"),
    make_option(c("--sumstats.frqID"), type="character", default="FRQ",
                help="Freq column name", metavar="character"),
    make_option(c("--strand.check"), type="numeric", default=0,
                help="Keep only non-ambiguous SNPs", metavar="numeric"),
    make_option(c("--param.file"), type="character", default=NULL,
                help="File of lambda and S parameters", metavar="character"),
    make_option(c("--ranking"), type="character", default="p-value",
                help="Method to rank loci", metavar="character"),

    make_option(c("--ld.shrink"), type="numeric",
                help="Indicator for shrinking LD matrix", metavar="numeric", default=0),
    make_option(c("--recomb.file"), type="character", default=NULL,
                help="File of genome-wide recombination rates", metavar="character"),
    make_option(c("--Ne"), type="numeric",
                help="Effective population size", metavar="numeric", default=0),

    make_option(c("--by.chr"), type="numeric", default=1,
                help="Logical indicating if bed files are split by chr",
                metavar="character"),
    make_option(c("--by.chr.sumstats"), type="character", default=0,
                help="Logical indicating if sumstats files are split by chr",
                metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print(opt)

logfile <- paste0(opt$beta.stem,".log")
tmp <- t(data.frame(opt))
rownames(tmp) <- names(opt)
write.table(tmp,file=logfile,quote=FALSE,col.names=FALSE)

w.prior <- as.numeric(strsplit( opt$w.prior, ',' )[[1]])

if( opt$by.chr.sumstats==0 ){
    sumstats <- fread( opt$sumstats, data.table=FALSE )
    sumstats.n <- sumstats[,opt$sumstats.nID]
    sumstats.se <- sumstats[,opt$sumstats.seID]
    sumstats.frq <- sumstats[,opt$sumstats.frqID]
    sigma2 <- median( 2*sumstats.n * sumstats.se * sumstats.se *
                      sumstats.frq * (1-sumstats.frq), na.rm=TRUE )
    snp.ptr <- which( colnames(sumstats)==opt$sumstats.snpID )
    allele1.ptr <- which( colnames(sumstats)==opt$sumstats.allele1ID )
    allele0.ptr <- which( colnames(sumstats)==opt$sumstats.allele0ID )
    beta.ptr <- which( colnames(sumstats)==opt$sumstats.betaID )
    sumstats <- sumstats[,c( snp.ptr, allele1.ptr, allele0.ptr, beta.ptr)]
    colnames(sumstats) <- c('SNP','ALLELE1','ALLELE0','BETA')
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
ld.ids <- as.character(read.table(opt$ld.ids)[,1])

if( opt$by.chr==0 ){
    ptr.bed <- BEDMatrix( opt$bfile, simple_names=TRUE )
    bim <- fread( paste(opt$bfile,'.bim',sep='' ) )
    ld.ids <- intersect( ld.ids, attributes(ptr.bed)[[3]][[1]] )
}

params <- read.table( opt$param.file, header=TRUE)
lambda.prior <- params$lambda.opt
S.prior <- params$S.opt

tmp <- strsplit( opt$prior, '/' )[[1]]
path <- paste(head(tmp,n=-1L),'/',sep='',collapse='')
lambda.ext <- paste( path, '/lambda/', sep='' )

for( chr in 1:22 ){
    if( opt$by.chr.sumstats!=0 ){
        sumstats <- fread( paste(opt$sumstats,chr,opt$by.chr.sumstats,sep=''),
                          data.table=FALSE )
        sumstats.n <- sumstats[,opt$sumstats.nID]
        if( chr==1 ){
            sumstats.se <- sumstats[,opt$sumstats.seID]
            sumstats.frq <- sumstats[,opt$sumstats.frqID]
            sigma2 <- median( 2*sumstats.n * sumstats.se * sumstats.se *
                      sumstats.frq * (1-sumstats.frq), na.rm=TRUE )
            snp.ptr <- which( colnames(sumstats)==opt$sumstats.snpID )
            allele1.ptr <- which( colnames(sumstats)==opt$sumstats.allele1ID )
            allele0.ptr <- which( colnames(sumstats)==opt$sumstats.allele0ID )
            beta.ptr <- which( colnames(sumstats)==opt$sumstats.betaID )
        }
        sumstats <- sumstats[,c( snp.ptr, allele1.ptr, allele0.ptr, beta.ptr)]
        colnames(sumstats) <- c('SNP','ALLELE1','ALLELE0','BETA')
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
    beta.prior <- fread(paste0(opt$prior,"_beta_bar_chr",chr,".txt.gz"))
    clump.ids <- unique(beta.prior$clump.id)
    beta.prior <- beta.list(beta.prior)
#    lambda.prior <- dget(paste(opt$prior,"_lambda_chr",chr,".Robj", sep=""))

    if( opt$ld.shrink==1 ){
        infile <- paste0(opt$recomb.file,'-',chr,'-final.txt.gz')
        recomb <- fread(infile)
    }
    if( opt$ld.shrink==0 ){
        recomb <- NULL
    }

    fits <- mclapply( 1:length(beta.prior),
                     function(i){
                         read.NonCentralFit.clump( beta.prior=beta.prior[[i]],
                                                  do.ld.shrink=opt$ld.shrink,
                                                  recomb=recomb, Ne=opt$Ne,
                                                  lambda.ext=lambda.ext,
                                                  w.prior=w.prior, sumstats=sumstats,
                                                  X.bed=ptr.bed, bim=bim,
                                                  ld.ids=ld.ids, by.chr=0,
                                                  sumstats.n=median(sumstats.n),
                                                  lambda.prior0=lambda.prior,
                                                  S.prior0=S.prior, ranking=opt$ranking,
                                                  strand.check=opt$strand.check )},
                     mc.cores=opt$n.cores )
#for( i in 1:length(beta.prior) ){
#    print(i)
#    tmp <- read.NonCentralFit.clump( beta.prior=beta.prior[[i]], do.ld.shrink=opt$ld.shrink, recomb=recomb, Ne=opt$Ne, lambda.ext=lambda.ext, w.prior=w.prior, sumstats=sumstats, X.bed=ptr.bed, bim=bim, ld.ids=ld.ids, by.chr=0, sumstats.n=median(sumstats.n), lambda.prior0=lambda.prior, S.prior0=S.prior, ranking=opt$ranking, strand.check=opt$strand.check )
#}

    names(fits) <- clump.ids
    fits.use <- which( !sapply( fits, is.null ) )
    fits <- fits[fits.use]
    beta.bar <- lapply( fits, getElement, 1 )
#    dput( beta.bar, paste(opt$beta.stem,"_beta_bar_chr",chr,".Robj", sep="") )

    outfile <- paste0(opt$beta.stem,"_beta_bar_chr",chr,".txt.gz")
    fwrite( beta.bar[[1]], outfile )
    for( i in 2:length(beta.bar) ){
        fwrite( beta.bar[[i]], outfile, append=TRUE )
    }
#    if( opt$precision==1 ){
#        lambda.fit <- lapply( fits, getElement, 2 )
#        dput( lambda.fit, paste(opt$beta.stem,"_lambda_chr",chr,".Robj", sep="" ) )
#    }
    kl.out <- t(simplify2array(lapply( fits, getElement, 2 )))
    kl.out <- data.frame( names(fits), kl.out )
    colnames(kl.out)[1] <- 'clumnp.id'
    outfile <- paste(opt$beta.stem,"_KLdist_chr",chr,".txt.gz", sep="")
    fwrite( kl.out, outfile )
}

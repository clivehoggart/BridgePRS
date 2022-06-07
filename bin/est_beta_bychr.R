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
    make_option(c("-s", "--sumstats"), type="character",
                help="Sum stats file", metavar="character"),
    make_option(c("--sumstats2"), type="character",
                help="Sum stats file pop2", metavar="character", default=NULL),
    make_option(c("--thinned.snplist"), type="character",
                help="List of SNPs to thin to", metavar="character", default=0),
    make_option(c("--n.max.locus"), type="numeric",
                help="Large loci are thinned", metavar="numeric", default=0),
    make_option(c("--ld.ids"), type="character",
                help="File listing IDs to use for estimating LD", metavar="character"),
    make_option(c("--param.file"), type="character", default=NULL,
                help="File of lambda and S parameters", metavar="character"),
    make_option(c("-b", "--beta.stem"), type="character", default=NULL,
                help="Out file stem", metavar="character"),
    make_option(c("--lambda","-l"), type="character", default="1",
                help="Comma delimited lambda values", metavar="character"),
    make_option(c("--S"), type="character", default="1",
                help="Comma delimited S values", metavar="character"),
    make_option(c("--precision"), type="logical", default=FALSE,
                help="Option to calculate precision matrix of each clump", metavar="character"),
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
    make_option(c("--ld.shrink"), type="numeric",
                help="Indicator for shrinking LD matrix", metavar="numeric", default=0),
    make_option(c("--recomb.file"), type="character", default=NULL,
                help="File of genome-wide recombination rates", metavar="character"),
    make_option(c("--Ne"), type="numeric",
                help="Effective population size", metavar="numeric", default=0),
    make_option(c("--by.chr"), type="logical", default=TRUE,
                help="Logical indicating if bed files are split by chr",
                metavar="character"),
    make_option(c("--by.chr.sumstats"), type="character", default=0,
                help="Logical indicating if sumstats files are split by chr",
                metavar="character"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

print(opt)

logfile <- paste0(opt$beta.stem,".log")
tmp <- t(data.frame(opt))
rownames(tmp) <- names(opt)
write.table(tmp,file=logfile,quote=FALSE,col.names=FALSE)

opt$beta.stem <- ifelse( is.null(opt$beta.stem), opt$clump.stem, opt$beta.stem )

if( opt$by.chr.sumstats==0 ){
    sumstats <- fread( opt$sumstats, data.table=FALSE )
    snp.ptr <- which( colnames(sumstats)==opt$sumstats.snpID )
    allele1.ptr <- which( colnames(sumstats)==opt$sumstats.allele1ID )
    allele0.ptr <- which( colnames(sumstats)==opt$sumstats.allele0ID )
    beta.ptr <- which( colnames(sumstats)==opt$sumstats.betaID )

    sumstats.n <- sumstats[,opt$sumstats.nID]
    sumstats.se <- sumstats[,opt$sumstats.seID]
    sumstats.frq <- sumstats[,opt$sumstats.frqID]
    sigma2 <- median( 2*sumstats.n * sumstats.se * sumstats.se *
                      sumstats.frq * (1-sumstats.frq), na.rm=TRUE )
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

lambda <- as.numeric(strsplit( opt$lambda, ',' )[[1]])
S <- as.numeric(strsplit( opt$S, ',' )[[1]])
precision <- opt$precision
if( !is.null(opt$param.file) ){
    params <- read.table( opt$param.file, header=TRUE)
    lambda <- params$lambda.opt
    S <- params$S.opt
}

ld.ids <- as.character(read.table(opt$ld.ids)[,1])

if( opt$by.chr==0 ){
    ptr.bed <- BEDMatrix( opt$bfile, simple_names=TRUE )
    bim <- fread( paste(opt$bfile,'.bim',sep='' ) )
    ld.ids <- intersect( ld.ids, attributes(ptr.bed)[[3]][[1]] )
}

for( chr in 1:22 ){
    if( opt$by.chr.sumstats!=0 ){
        sumstats <- fread( paste(opt$sumstats,chr,opt$by.chr.sumstats,sep=''), data.table=FALSE )
        if( chr==1 ){
            sumstats.n <- sumstats[,opt$sumstats.nID]
            sumstats.se <- sumstats[,opt$sumstats.seID]
            sumstats.frq <- sumstats[,opt$sumstats.frqID]
            sigma2 <- median( 2*sumstats.n * sumstats.se * sumstats.se *
                              sumstats.frq * (1-sumstats.frq), na.rm=TRUE )
        }
        snp.ptr <- which( colnames(sumstats)==opt$sumstats.snpID )
        allele1.ptr <- which( colnames(sumstats)==opt$sumstats.allele1ID )
        allele0.ptr <- which( colnames(sumstats)==opt$sumstats.allele0ID )
        beta.ptr <- which( colnames(sumstats)==opt$sumstats.betaID )
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
    infile <- paste(opt$clump.stem,'_',chr,'.clumped.gz',sep='')
    clump <- read.table(infile,header=TRUE,stringsAsFactors=FALSE)
    clump.use <- which( clump$P <= 1 )

    if( opt$ld.shrink==1 ){
        infile <- paste0(opt$recomb.file,'-',chr,'-final.txt.gz')
        recomb <- fread(infile)
    }
    if( opt$ld.shrink==0 ){
        recomb <- NULL
    }

    if( opt$thinned.snplist!=0 ){
        thinned.snplist <- fread(opt$thinned.snplist)$V2
        if( !is.null(opt$sumstats2) ){
            infile <- paste0( opt$sumstats2, chr, opt$by.chr.sumstats )
            sumstats2 <- fread(infile)
            snp.ptr <- which( colnames(sumstats2)==opt$sumstats.snpID )
            allele1.ptr <- which( colnames(sumstats2)==opt$sumstats.allele1ID )
            allele0.ptr <- which( colnames(sumstats2)==opt$sumstats.allele0ID )
            beta.ptr <- which( colnames(sumstats2)==opt$sumstats.betaID )
            sumstats2 <- sumstats2[,c( snp.ptr, allele1.ptr, allele0.ptr, beta.ptr)]
            colnames(sumstats2) <- c('SNP','ALLELE1','ALLELE0','BETA')
            clump <- thin.big.loci( clump, thinned.snplist, opt$n.max.locus, sumstats2 )
        }else{
            clump <- thin.big.loci( clump, thinned.snplist, opt$n.max.locus )
        }
    }

    tmp <- strsplit( opt$beta.stem, '/' )[[1]]
    path <- paste(head(tmp,n=-1L),'/',sep='',collapse='')
    fits <- mclapply( 1:length(clump.use),
                     function(i){
                         read.fit.clump( clump.i=clump[clump.use[i],],
                                        do.ld.shrink=opt$ld.shrink,
                                        sumstats=sumstats,
                                        recomb=recomb, Ne=opt$Ne,
                                        X.bed=ptr.bed, bim=bim, ld.ids=ld.ids,
                                        S=S, l=lambda, precision=precision, by.chr=0,
                                        beta.stem=path,
                                        strand.check=opt$strand.check )},
                     mc.cores=as.numeric(opt$n.cores) )

#    for( i in 1:length(clump.use) ){
#        print(i)
#        tmp <- read.fit.clump( clump.i=clump[clump.use[i],],
#                               do.ld.shrink=opt$ld.shrink, Ne=opt$Ne,
#                               sumstats=sumstats, bim=bim, ld.ids=ld.ids, S=S,
#                               recomb=recomb, precision=precision,
#                               X.bed=ptr.bed, l=lambda,
#                               by.chr=0, beta.stem=path, strand.check=opt$strand.check )
#    }
#fits <- read.fit.clump( clump.i=clump[clump.use[i],], do.ld.shrink=opt$ld.shrink, sumstats=sumstats, recomb=recomb, Ne=opt$Ne, X.bed=ptr.bed, bim=bim, ld.ids=ld.ids, S=S, l=lambda, precision=precision, by.chr=0, beta.stem=path )
    beta.bar <- lapply( fits, getElement, 1 )
    names(beta.bar) <- clump$SNP[clump.use]

    ptr.qc <- which(sapply(beta.bar,length)!=0)
    beta.bar <- beta.bar[ptr.qc]
#    dput( beta.bar, paste(opt$beta.stem,"_beta_bar_chr",chr,".Robj", sep="") )

    outfile <- paste(opt$beta.stem,"_beta_bar_chr",chr,".txt.gz", sep="")
    fwrite( beta.bar[[1]], outfile )
    for( i in 2:length(beta.bar) ){
        fwrite( beta.bar[[i]], outfile, append=TRUE )
    }

#    kl.out <- t(simplify2array(lapply( fits, getElement, 2 )))
#    if( !is.na(kl.out[1]) ){
#        kl.out <- data.frame( clump$SNP[clump.use], kl.out )
#        colnames(kl.out)[1] <- 'clumnp.id'
#        outfile <- paste(opt$beta.stem,"_KLdist_chr",chr,".txt.gz", sep="")
#        fwrite( kl.out, outfile )
#    }
}

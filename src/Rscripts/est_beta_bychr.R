library(BEDMatrix)
library(MASS)
library(parallel)
library(data.table)
library("optparse")
options(stringsAsFactors=FALSE)


# toupper

option_list = list(
    make_option(c("--fpath"), type="character", default=NULL,help="Function File Path", metavar="character"),
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
                help="Option to calculate precision matrix of each clump",
                metavar="character"),
    make_option(c("--n.cores"), type="numeric", default=1,
                help="Number of processors for mclapply to use", metavar="character"),
    make_option(c("--bfile"), type="character",
                help="Plink file to estimate LD", metavar="character"),
    make_option(c("--sumstats.snpID"), type="character", default="SNP",
                help="SNP column name", metavar="character"),
    make_option(c("--sumstats.betaID"), type="character", default="BETA",
                help="Beta column name", metavar="character"),
    make_option(c("--sumstats.allele0ID"), type="character", default="ALLELE0",
                help="Allele 0 column name", metavar="character"),
    make_option(c("--sumstats.allele1ID"), type="character", default="ALLELE1",
                help="Allele 1 column name", metavar="character"),
    make_option(c("--strand.check"), type="numeric", default=0,
                help="Keep only non-ambiguous SNPs", metavar="numeric"),
    make_option(c("--by.chr"), type="numeric", default=1,
                help="Logical indicating if bed files are split by chr",
                metavar="character"),
    make_option(c("--by.chr.sumstats"), type="character", default=0,help="Logical indicating if sumstats files are split by chr",metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

print(opt)

source(opt$fpath)


logfile <- paste0(opt$beta.stem,".log")
tmp <- t(data.frame(opt))

rownames(tmp) <- names(opt)
write.table(tmp,file=logfile,quote=FALSE,col.names=FALSE)


opt$beta.stem <- ifelse( is.null(opt$beta.stem), opt$clump.stem, opt$beta.stem )

lambda <- as.numeric(strsplit( opt$lambda, ',' )[[1]])
S <- as.numeric(strsplit( opt$S, ',' )[[1]])
precision <- opt$precision

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
    if( !is.null(opt$param.file) ){
        params <- read.table( opt$param.file, header=TRUE)
        lambda <- params$lambda.opt
        p <- params$p.opt
        S <- params$S.opt
#        out <- data.frame(S,lambda,p)
#        colnames(out) <- c("S.opt","lambda.opt","p.opt")
#        write.table( out, opt$param.file, row.names=FALSE, quote=FALSE )
    }
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
        if( chr==1 ){
            if( !is.null(opt$param.file) ){
                params <- read.table( opt$param.file, header=TRUE)
                lambda <- params$lambda.opt
                p <- params$p.opt
                S <- params$S.opt
#                out <- data.frame(S,lambda,p)
#                colnames(out) <- c("S.opt","lambda.opt","p.opt")
#                write.table( out, opt$param.file, row.names=FALSE, quote=FALSE )
            }
        }

        snp.ptr <- which( colnames(sumstats)==opt$sumstats.snpID )
        allele1.ptr <- which( colnames(sumstats)==opt$sumstats.allele1ID )
        allele0.ptr <- which( colnames(sumstats)==opt$sumstats.allele0ID )
        beta.ptr <- which( colnames(sumstats)==opt$sumstats.betaID )
        #print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXi44')
        #print(snp.ptr)
        #print(allele1.ptr)
        #print(allele0.ptr)
        #print(beta.ptr)

        sumstats <- sumstats[,c( snp.ptr, allele1.ptr, allele0.ptr, beta.ptr)]

        colnames(sumstats) <- c('SNP','ALLELE1','ALLELE0','BETA')
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
    infile <- paste(opt$clump.stem,'_',chr,'.clumped.gz',sep='')
    clump <- read.table(infile,header=TRUE,stringsAsFactors=FALSE)
    clump.use <- which( clump$P <= 1 )

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
                                        sumstats=sumstats, ld.ids=ld.ids,
                                        X.bed=ptr.bed, bim=bim,
                                        l=lambda, S=S, precision=precision, by.chr=0,
                                        beta.stem=path,
                                        strand.check=opt$strand.check )},
                     mc.cores=as.numeric(opt$n.cores) )

#    for( i in 1:length(clump.use) ){
#        print(i)
#        tmp <- read.fit.clump( clump.i=clump[clump.use[i],],
#                                        sumstats=sumstats,
#                                        X.bed=ptr.bed, bim=bim, ld.ids=ld.ids,
#                                        S=S, l=lambda, precision=precision, by.chr=0,
#                                        beta.stem=path,
#                                        strand.check=opt$strand.check )
#    }
    beta.bar <- lapply( fits, getElement, 1 )

    names(beta.bar) <- clump$SNP[clump.use]

    ptr.qc <- which(sapply(beta.bar,length)!=0)
    beta.bar <- beta.bar[ptr.qc]

    outfile <- paste(opt$beta.stem,"_beta_bar_chr",chr,".txt.gz", sep="")
    fwrite(rbindlist(beta.bar), outfile)
 #   fwrite( beta.bar[[1]], outfile )
 #   for( i in 2:length(beta.bar) ){
 #       fwrite( beta.bar[[i]], outfile, append=TRUE )
 #   }

#    kl.out <- t(simplify2array(lapply( fits, getElement, 2 )))
#    if( !is.na(kl.out[1]) ){
#        kl.out <- data.frame( clump$SNP[clump.use], kl.out )
#        colnames(kl.out)[1] <- 'clumnp.id'
#        outfile <- paste(opt$beta.stem,"_KLdist_chr",chr,".txt.gz", sep="")
#        fwrite( kl.out, outfile )
#    }
}
if( !is.null(warnings()) ){
    print(warnings())
}

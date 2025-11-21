library(data.table)
library("optparse")
options(stringsAsFactors=FALSE)
library(BEDMatrix)

blocks.stats <- function( block.i, sumstats, models ){
    Sigma.prs <- list()
    betatXtY.2 <- list()
        betatXtY.3 <- list()
    block.snps.all <-  unlist(strsplit( block.i$SNPS, '\\|' ))
    block.snps1 <- vector()
    for( k in 1:n.models ){
        block.snps.k <- intersect( rownames(models[[k]]), block.snps.all )
        if( length(block.snps.k)>0 ){
            ptr.prs <- match( block.snps.k, rownames(models[[k]]) )
            ptr.ss <- match( block.snps.k, sumstats$SNP )
            ref.stats <- est.ref.stats( block.snps.k, ld.ids, X.bed, bim,
                                       sumstats$ALLELE1[ptr.ss],
                                       sumstats$ALLELE0[ptr.ss],
                                       opt$strand.check, n.eff=FALSE )
            beta <- models[[k]][ptr.prs,,drop=FALSE]
            Sigma.prs[[k]] <- t(beta) %*% ref.stats$ld %*% beta
            betatXtY.2[[k]] <- t(beta) %*% as.matrix(sumstats$XtY.2[ptr.ss])
            betatXtY.3[[k]] <- t(beta) %*% as.matrix(sumstats$XtY.3[ptr.ss])
            block.snps1 <- union( block.snps1, block.snps.k )
        }
    }
    return( list( block.snps1, Sigma.prs, betatXtY.2, betatXtY.3 ) )
}

# toupper

option_list = list(
    make_option(c("--fpath"), type="character", default=NULL,
                help="Function File Path", metavar="character"),
    make_option("--stage1", type="character",
                help="Stage 1 betas", metavar="character"),
    make_option("--stage2", type="character", default=NULL,
                help="Stage 2 betas", metavar="character"),
    make_option(c("--p.thresh"), type="character",
                default="1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8",
                help="P-value threshold to select clumps",
                metavar="character"),
    make_option(c("--workdir"), type="character",
                help="Top level working dir", metavar="character"),
    make_option(c("--blockdir"), type="character",
                help="Block dir", metavar="character"),
    make_option(c("--bfile"), type="character",
                help="Plink file to estimate LD", metavar="character"),
    make_option(c("--ld.ids"), type="character",
                help="File listing IDs to use for estimating LD", metavar="character"),
    make_option(c("--by.chr"), type="numeric", default=1,
                help="Logical indicating if bed files are split by chr",
                metavar="character"),
    make_option(c("--strand.check"), type="numeric", default=0,
                help="Keep only non-ambiguous SNPs", metavar="numeric"),
    make_option(c("--fold"), type="numeric", default=0,
                help="Fold", metavar="numeric"),
    make_option(c("--n.cores"), type="numeric", default=1,
                help="Number of processors for mclapply to use", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print(opt)

source(opt$fpath)

p.thresh <- as.numeric(strsplit( opt$p.thresh, ',' )[[1]])
ld.ids <- as.character(read.table(opt$ld.ids)[,2])
if( opt$by.chr==0 ){
    ptr.bed <- BEDMatrix( opt$bfile, simple_names=TRUE )
    bim <- fread( paste(opt$bfile,'.bim',sep='' ) )
    ld.ids <- intersect( ld.ids, attributes(ptr.bed)[[3]][[1]] )
}

af <- vector()
models <- list()
for( chr in 1:22 ){
    if( opt$by.chr==1 ){
        X.bed <- BEDMatrix( paste(opt$bfile,chr,sep=''), simple_names=TRUE )
        bim <- fread( paste(opt$bfile,chr,'.bim',sep='' ) )
        ld.ids <- intersect( ld.ids, attributes(X.bed)[[3]][[1]] )
    }

    blockfile <- paste(opt$blockdir,'/chr',chr,'.blocks.det.gz',sep='')
    blocks <- fread( blockfile, data.table=FALSE )
    sumstatfile <- paste(opt$workdir,'/fold',opt$fold,
                         '/sumstat_subset/chr',chr,'.dat.gz',sep='')
    sumstats <- fread( sumstatfile, data.table=FALSE )

    infile <- paste( opt$stage1,"_beta_bar_chr",chr,".txt.gz", sep="")
    stage1 <- fread( infile, data.table=FALSE )
    if( !is.null(opt$stage2) ){
        infile <- paste(opt$stage2,"_beta_bar_chr",chr,".txt.gz", sep="")
        stage2 <- fread( infile, data.table=FALSE )
        infile <- paste0(opt$stage2,"_KLdist_chr",chr,".txt.gz")
#        kl.metric <- fread(infile, data.table=FALSE)[,-1]
        snps <- union( stage1$snp, stage2$snp )
    }else{
        snps <- stage1$snp
    }
    ptr <- match( snps, sumstats$SNP )
    chr.alleles <- data.frame( sumstats$ALLELE1[ptr], sumstats$ALLELE0[ptr] )
    colnames(chr.alleles) <- c('ALLELE1','ALLELE0')

    ptr.betas <- grep('beta.bar',colnames(stage1))
    tmp <- strsplit(colnames(stage1)[ptr.betas],'beta.bar_')
    betas <- sapply(tmp,getElement,2)
    models <- as.data.frame(matrix(ncol=length(p.thresh)*length(betas),
                                    nrow=length(snps),data=0))
    cnames <- vector()
    for( i in 1:length(p.thresh) ){
        cnames <- c( cnames, paste('beta',betas,p.thresh[i],sep="_") )
    }
    colnames(models) <- cnames
    ptr.ss <- match( stage1$snp, sumstats$SNP )
    swtch <- allele.check( sumstats$ALLELE1[ptr.ss], sumstats$ALLELE0[ptr.ss],
                          stage1$effect.allele, stage1$ref.allele,
                          opt$strand.check )
    for( i in 1:length(p.thresh) ){
        ptr.stage1 <- which( stage1$p.value <  p.thresh[i] )
        ptr.row <- match( stage1$snp[ptr.stage1], snps )
        ptr.col <- (1:length(betas)) + (i-1)*length(betas)
        models[ptr.row,ptr.col] <- stage1[ptr.stage1,ptr.betas] * swtch[ptr.stage1]
    }
    rownames(models) <- snps
    models <- as.matrix(models)

    if( !is.null(opt$stage2) ){
        ptr.betas <- grep('beta.bar',colnames(stage2))
        tmp <- strsplit(colnames(stage2)[ptr.betas],'beta.bar_')
        betas <- sapply(tmp,getElement,2)
        models2 <- as.data.frame(matrix(ncol=length(p.thresh)*length(betas),
                                        nrow=length(snps),data=0))
        cnames <- vector()
        for( i in 1:length(p.thresh) ){
            cnames <- c( cnames, paste('beta',betas,p.thresh[i],sep="_") )
        }
        colnames(models2) <- cnames
        ptr.ss <- match( stage2$snp, sumstats$SNP )
        swtch <- allele.check( sumstats$ALLELE1[ptr.ss], sumstats$ALLELE0[ptr.ss],
                              stage2$effect.allele, stage2$ref.allele,
                              opt$strand.check )
        for( i in 1:length(p.thresh) ){
            ptr.stage2 <- which( stage2$p.value <  p.thresh[i] )
            ptr.row <- match( stage2$snp[ptr.stage2], snps )
            ptr.col <- (1:length(betas)) + (i-1)*length(betas)
            models2[ptr.row,ptr.col] <- stage2[ptr.stage2,ptr.betas] * swtch[ptr.stage2]
        }
        rownames(models2) <- snps
        models <- as.matrix(cbind(models,models2))
    }
    n.models <- length(models)
#    fwrite( all.models, outfile )
    if( chr==1 ){
        genome.alleles <- chr.alleles
        genome.models <- models
        kk <- ncol(models)
        Sigma.prs <- matrix( ncol=kk, nrow=kk, data=0 )
        betatXtY.2 <- matrix( ncol=1, nrow=kk, data=0 )
        betatXtY.3 <- matrix( ncol=1, nrow=kk, data=0 )
    }else{
        genome.alleles <- rbind( genome.alleles, chr.alleles )
        genome.models <- rbind( genome.models, models )
    }

    all.block.snps <- vector()
    k <- n.models
    for( i in 1:nrow(blocks) ){
        block.snps.all <-  unlist(strsplit( blocks$SNPS[i], '\\|' ))
#        block.snps <- list()
#        block.snps1 <- vector()
#        for( k in 1:n.models ){
        block.snps <- intersect( rownames(models), block.snps.all )
        if( length(block.snps) > 0 ){
            ptr.prs <- match( block.snps, rownames(models) )
            ptr.ss <- match( block.snps, sumstats$SNP )
            ref.stats <- est.ref.stats( block.snps, ld.ids, X.bed, bim,
                                       sumstats$ALLELE1[ptr.ss],
                                       sumstats$ALLELE0[ptr.ss],
                                       opt$strand.check, n.eff=FALSE )
            beta <- models[ptr.prs,,drop=FALSE]
            Sigma.prs <- Sigma.prs + t(beta) %*% ref.stats$ld %*% beta
            betatXtY.2 <- betatXtY.2 + t(beta) %*% as.matrix(sumstats$XtY.2[ptr.ss])
            betatXtY.3 <- betatXtY.3 + t(beta) %*% as.matrix(sumstats$XtY.3[ptr.ss])
            all.block.snps <- union( all.block.snps, block.snps )
#            }
#            all.block.snps <- c( all.block.snps, block.snps1 )
        }
    }
    single.snps <- setdiff( snps, all.block.snps )
    for( i in 1:length(single.snps) ){
        ptr.ss <- match( single.snps[i], sumstats$SNP )
#        for( k in 1:n.models ){
            ptr.prs <- match( single.snps[i], rownames(models) )
            ref.stats <- est.ref.stats( single.snps[i], ld.ids, X.bed, bim,
                                       sumstats$ALLELE1[ptr.ss],
                                       sumstats$ALLELE0[ptr.ss],
                                       opt$strand.check, n.eff=FALSE )
            beta <- models[ptr.prs,,drop=FALSE]
            Sigma.prs <- Sigma.prs +  t(beta) %*% ref.stats$ld %*% beta
            betatXtY.2 <- betatXtY.2 + t(beta) %*% as.matrix(sumstats$XtY.2[ptr.ss])
            betatXtY.3 <- betatXtY.3 + t(beta) %*% as.matrix(sumstats$XtY.3[ptr.ss])
#        }
    }
}

save.image(paste0( opt$workdir,'/fold',opt$fold,'/debug.RData'))
#    tmp <- mclapply( 1:nrow(blocks),
#                    function(i){
#                        block.stats( block.i=blocks[i,], sumstats, models )},
#                    mc.cores=as.numeric(opt$n.cores) )

#for( kk in 1:10 ){
#    load(paste0('/sc/arion/projects/psychgen/projects/prs/cross_population_prs_development/quick_ridge/results/sumstats/hm/50/EUR/fold',kk,'/debug.RData'))

R2.indiv <- 2*log(abs( betatXtY.2 + betatXtY.3 )) - log(diag(Sigma.prs))
s2 <- order( R2.indiv, decreasing=TRUE )
tmp <- strsplit( colnames(models)[s2[1]], '_' )
S.opt <- as.numeric(tmp[[1]][3])
lambda.opt <- as.numeric(tmp[[1]][2])
p.opt <- as.numeric(tmp[[1]][4])
write.table( data.frame(S.opt,lambda.opt,p.opt),
            paste0( opt$workdir,'/fold',opt$fold,'/prs_sumstats/best_model_params.dat' ),
            row.names=FALSE, quote=FALSE )
#}

write.table( betatXtY.2, paste0( opt$workdir,'/fold',opt$fold,'/prs_sumstats/betatXtY_2.dat' ),
            row.names=TRUE, col.names=FALSE, quote=FALSE )
write.table( betatXtY.3, paste0( opt$workdir,'/fold',opt$fold,'/prs_sumstats/betatXtY_3.dat' ),
            row.names=TRUE, col.names=FALSE, quote=FALSE )
fwrite( Sigma.prs, paste0( opt$workdir,'/fold',opt$fold,'/prs_sumstats/Sigma_prs.dat' ) )
fwrite( data.frame( rownames(genome.models), genome.alleles, genome.models ),
       paste0( opt$workdir,'/fold',opt$fold,'/prs_sumstats/snp_weights.dat' ) )
#}
#save.image(paste0( opt$workdir,'/fold',opt$fold,'/debug.RData'))

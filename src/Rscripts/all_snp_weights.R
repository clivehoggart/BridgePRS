library(data.table)
library("optparse")
options(stringsAsFactors=FALSE)
library(BEDMatrix)

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
    make_option(c("--bfile"), type="character",
                help="Plink file to estimate LD", metavar="character"),
    make_option(c("--ld.ids"), type="character",
                help="File listing IDs to use for estimating LD", metavar="character"),
    make_option(c("--by.chr"), type="numeric", default=1,
                help="Logical indicating if bed files are split by chr",
                metavar="character"),
    make_option(c("--strand.check"), type="numeric", default=0,
                help="Keep only non-ambiguous SNPs", metavar="numeric"),
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

for( chr in 1:22 ){
    if( opt$by.chr==1 ){
        X.bed <- BEDMatrix( paste(opt$bfile,chr,sep=''), simple_names=TRUE )
        bim <- fread( paste(opt$bfile,chr,'.bim',sep='' ) )
        ld.ids <- intersect( ld.ids, attributes(X.bed)[[3]][[1]] )
    }

    blockfile <- paste(opt$workdir,'/blocks/chr',chr,'.blocks.det.gz',sep='')
    blocks <- fread( blockfile, data.table=FALSE )
    sumstatfile <- paste(opt$workdir,'/sumstat_subset/chr',chr,'.dat.gz',sep='')
    sumstats <- fread( sumstatfile, data.table=FALSE )

    infile <- paste( opt$stage1,"_beta_bar_chr",chr,".txt.gz", sep="")
    stage1 <- fread( infile, data.table=FALSE )
    if( !is.null(opt$stage2) ){
        infile <- paste(opt$stage2,"_beta_bar_chr",chr,".txt.gz", sep="")
        stage2 <- fread( infile, data.table=FALSE )
        infile <- paste0(opt$stage2,"_KLdist_chr",chr,".txt.gz")
        kl.metric <- fread(infile, data.table=FALSE)[,-1]
        snps <- union( stage1$snp, stage2$snp )
    }else{
        snps <- stage1$snp
    }

    ptr.betas <- grep('beta.bar',colnames(stage1))
    tmp <- strsplit(colnames(stage1)[ptr.betas],'beta.bar_')
    betas <- sapply(tmp,getElement,2)

    ptr <- match( snps, sumstats$SNP )
    chr.alleles <- data.frame( sumstats$ALLELE1[ptr], sumstats$ALLELE0[ptr] )
    colnames(chr.alleles) <- c('ALLELE1','ALLELE0')

    stage11 <- as.data.frame(matrix(ncol=length(p.thresh)*length(betas),
                                    nrow=length(snps),data=0))
    cnames <- vector()
    for( i in 1:length(p.thresh) ){
        cnames <- c( cnames, paste('beta',betas,p.thresh[i],sep="_") )
    }
    colnames(stage11) <- cnames
    ptr.ss <- match( stage1$snp, sumstats$SNP )
    swtch <- allele.check( sumstats$ALLELE1[ptr.ss], sumstats$ALLELE0[ptr.ss],
                          stage1$effect.allele, stage1$ref.allele,
                          opt$strand.check )
    for( i in 1:length(p.thresh) ){
        ptr.stage1 <- which( stage1$p.value <  p.thresh[i] )
        ptr.row <- match( stage1$snp[ptr.stage1], snps )
        ptr.col <- (1:length(betas)) + (i-1)*length(betas)
        stage11[ptr.row,ptr.col] <- stage1[ptr.stage1,ptr.betas] * swtch[ptr.stage1]
    }

    if( !is.null(opt$stage2) ){
        if( chr==1 ){
            clump.id <- unique(stage2$clump.id)
            p.clump <- stage2$p.value[match( clump.id, stage2$clump.id )]
            F <- ecdf(p.clump)
            q <- F(10^(-(1:8)))
            kl.thresh <- apply( kl.metric, 2, quantile, 1-q, na.rm=TRUE )
            n.thresh <- nrow(kl.thresh)
            nme.thresh <- rownames(kl.thresh)
        }
        ptr.betas <- grep('beta.bar',colnames(stage2))
        tmp <- strsplit(colnames(stage2)[ptr.betas],'beta.bar_')
        betas <- sapply(tmp,getElement,2)
        stage22 <- as.data.frame(matrix(ncol=length(kl.thresh)*length(betas),nrow=length(snps),data=0))
        cnames <- vector()
        for( i in 1:length(kl.thresh) ){
            cnames <- c( cnames, paste('beta',betas,kl.thresh[i],sep="_") )
        }
        colnames(stage22) <- cnames
        ptr.ss <- match( stage2$snp, sumstats$SNP )
        swtch <- allele.check( sumstats$ALLELE1[ptr.ss], sumstats$ALLELE0[ptr.ss],
                              stage2$effect.allele, stage2$ref.allele,
                              opt$strand.check )
        for( i in 1:length(kl.thresh) ){
            ptr.stage2 <- which( stage2$p.value <  kl.thresh[i] )
            ptr.row <- match( stage2$snp[ptr.stage2], snps )
            ptr.col <- (1:length(betas)) + (i-1)*length(betas)
            stage22[ptr.row,ptr.col] <- stage2[ptr.stage2,ptr.betas] * swtch[ptr.stage2]
        }
        all.models <- as.matrix(cbind(stage11,stage22))
    }else{
        all.models <- as.matrix(stage11)
    }
    rownames(all.models) <- snps
    n.models <- ncol(all.models)
#    fwrite( all.models, outfile )
    if( chr==1 ){
        genome.alleles <- chr.alleles
        genome.all.models <- all.models
        Sigma.prs <- matrix( ncol=n.models, nrow=n.models, data=0 )
        betatXtY.2 <- matrix( ncol=1, nrow=n.models, data=0 )
        betatXtY.3 <- matrix( ncol=1, nrow=n.models, data=0 )
    }else{
        genome.alleles <- rbind( genome.alleles, chr.alleles )
        genome.all.models <- rbind( genome.all.models, all.models )
    }
    genome.all.models <- standardise(genome.all.models)

    all.block.snps <- vector()
    for( i in 1:nrow(blocks) ){
        block.snps <-  unlist(strsplit( blocks$SNPS[i], '\\|' ))
        block.snps <- intersect( rownames(all.models), block.snps )
#        print(c(i,length(block.snps)))}
        if( length(block.snps)>0 ){
            all.block.snps <- c( all.block.snps, block.snps )
            ptr.ss <- match( block.snps, sumstats$SNP )
            ptr.prs <- match( block.snps, rownames(all.models) )
            ref.stats <- est.ref.stats( block.snps, ld.ids, X.bed, bim,
                                       sumstats$ALLELE1[ptr.ss],
                                       sumstats$ALLELE0[ptr.ss],
                                       opt$strand.check, n.eff=FALSE )
            beta <- all.models[ptr.prs,,drop=FALSE]
            Sigma.prs <- Sigma.prs + t(beta) %*% ref.stats$ld %*% beta
            betatXtY.2 <- betatXtY.2 + t(beta) %*% as.matrix(sumstats$XtY.2[ptr.ss])
            betatXtY.3 <- betatXtY.3 + t(beta) %*% as.matrix(sumstats$XtY.3[ptr.ss])
        }
    }
    single.snps <- setdiff( rownames(all.models), all.block.snps )
    for( i in 1:length(single.snps) ){
        ptr.ss <- match( single.snps[i], sumstats$SNP )
        ptr.prs <- match( single.snps[i], rownames(all.models) )
        ref.stats <- est.ref.stats( single.snps[i], ld.ids, X.bed, bim,
                                   sumstats$ALLELE1[ptr.ss], sumstats$ALLELE0[ptr.ss],
                                   opt$strand.check, n.eff=FALSE )
        beta <- all.models[ptr.prs,,drop=FALSE]
        Sigma.prs <- Sigma.prs +  t(beta) %*% ref.stats$ld %*% beta
        betatXtY.2 <- betatXtY.2 + t(beta) %*% as.matrix(sumstats$XtY.2[ptr.ss])
        betatXtY.3 <- betatXtY.3 + t(beta) %*% as.matrix(sumstats$XtY.3[ptr.ss])
    }
}

lambda <- c((1:9)*1e-6,(1:9)*1e-5,(1:9)*1e-4,(1:9)*1e-3,(1:9)*1e-2,(1:9)*1e-1,1:9)
prs.weights <- matrix( ncol=length(lambda), nrow=n.models )
R2.ensembl <- vector()
for( k in 1:length(lambda) ){
    prs.weights[,k] <- solve( diag(lambda[k],n.models) + Sigma.prs ) %*% betatXtY.2
#    R2.ensembl[k] <- (t(prs.weights[,k,drop=FALSE]) %*% betatXtY.3)^2 / t(prs.weights[,k,drop=FALSE]) %*% Sigma.prs %*% prs.weights[,k,drop=FALSE]
}
R2.ensembl <- (t(prs.weights) %*% betatXtY.3)^2 / diag(t(prs.weights) %*% Sigma.prs %*% prs.weights)
s1 <- order( R2.ensembl, decreasing=TRUE )

ensembl.model <- genome.all.models %*% prs.weights[,s1[1]]

write.table( data.frame( genome.alleles, ensembl.model ),
            paste( opt$workdir, '_snp_weights_weighted_model.dat', sep='' ),
            col.names=FALSE, quote=FALSE )

#cbind(lambda,R2.ensembl/max(R2.ensembl))
#s[1:4]

R2.indiv <- (betatXtY.2 + betatXtY.3)^2 / diag(Sigma.prs)
s2 <- order( R2.indiv, decreasing=TRUE )
colnames(all.models)[s2[1:10]]

tmp <- strsplit( colnames(all.models)[s2 [1]], '_' )
S.opt <- as.numeric(tmp[[1]][3])
lambda.opt <- as.numeric(tmp[[1]][2])
p.opt <- as.numeric(tmp[[1]][4])
write.table( data.frame(S.opt,lambda.opt,p.opt),
            paste( opt$workdir, '_best_model_params.dat', sep='' ),
            row.names=FALSE, quote=FALSE )
write.table( data.frame( genome.alleles, genome.all.models[,s2[1]] ),
            paste( opt$workdir, '_snp_weights_best_model.dat', sep='' ),
            col.names=FALSE, quote=FALSE )

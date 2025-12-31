library(data.table)
library("optparse")

enet.R2 <- function( lambda, n, Sigma.prs, betatXtY.2, betatXtY.3 ){
    n.prs <- nrow(Sigma.prs)
    prs.weights <- solve( diag(lambda,n.prs) + n*Sigma.prs ) %*% betatXtY.2
    log.R2.ensembl <- 2*log(abs(t(prs.weights) %*% betatXtY.3)) -
        log(t(prs.weights) %*% Sigma.prs %*% prs.weights)
    return(log.R2.ensembl)
}
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


option_list = list(
    make_option(c("--N.pop"), type="numeric", default=0,
                help="GWAS sample size", metavar="numeric"),
    make_option(c("--prop.train"), type="numeric", default=0,
                help="GWAS training sample proportion", metavar="numeric"),
    make_option(c("--prop.test"), type="numeric", default=0,
                help="Testing sample proportion", metavar="numeric"),
    make_option(c("--n.folds"), type="numeric", default=0,
                help="cv folds", metavar="numeric"),
    make_option("--workdir", type="character", default=NULL,
                help="file path", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print(opt)

#opt <- list()
#opt$N.pop <- 4800
#opt$prop.train <- 0.6
#opt$prop.test <- 0.3
#opt$n.folds <- 10
#opt$workdir <- "/sc/arion/projects/psychgen/projects/prs/cross_population_prs_development/quick_ridge/results/sumstats/hm/50/AFR"

n.test <- opt$N.pop * opt$prop.test
n.test.all <- opt$N.pop * (1 - opt$prop.train)
R2.model <- matrix( ncol=3, nrow=opt$n.folds )
ensembl.model <- list()
ptr.use <- list()
lambda <- matrix( ncol=3, nrow=opt$n.folds )
betatXtY.2 <- list(length=opt$n.folds)
betatXtY.3 <- list(length=opt$n.folds)
Sigma.prs <- list(length=opt$n.folds)
genome.models <- list(length=opt$n.folds)
alleles <- list(length=opt$n.folds)
model.names <- list(length=opt$n.folds)
snps <- vector()
R2.indiv <- list()
all.model.names <- vector()
for( kk in 1:opt$n.folds ){
    betatXtY.2[[kk]] <- list(length=3)
    betatXtY.3[[kk]] <- list(length=3)
    Sigma.prs[[kk]] <- list(length=3)
    genome.models[[kk]] <- list(length=3)
    model.names[[kk]] <- list(length=3)

    betatXtY.22 <- fread( paste0( opt$workdir,'/fold',kk,'/prs_sumstats/betatXtY_2.dat' ) )
    betatXtY.33 <- fread( paste0( opt$workdir,'/fold',kk,'/prs_sumstats/betatXtY_3.dat' ) )
    Sigma.prs1 <- fread( paste0( opt$workdir,'/fold',kk,'/prs_sumstats/Sigma_prs.dat' ) )
    genome.models1 <- fread( paste0( opt$workdir,'/fold',kk,'/prs_sumstats/snp_weights.dat' ) )

    alleles[[kk]] <- data.frame( genome.models1$rownames, genome.models1$ALLELE1, genome.models1$ALLELE0 )
    snps1 <- genome.models1$rownames
    snps <- union( snps, snps1 )
    genome.models1 <- as.matrix(genome.models1[,-(1:3)])

    model.type <- sapply(strsplit( betatXtY.22$V1, "_" ),length)
    ptr.models <- list(lengtgh=3)
    ptr.models[[1]] <- which( model.type==4 )
    ptr.models[[2]] <- which( model.type==3 )
    ptr.models[[3]] <-  1:nrow(betatXtY.22)
    if( length(ptr.models[[1]])>0 & length(ptr.models[[2]])>0 ){
        n.models <- 3
    }else{
        n.models <- 1
    }
    for( k in 1:n.models ){
        Sigma.prs[[kk]][[k]] <- as.matrix(Sigma.prs1)[ ptr.models[[k]], ptr.models[[k]] ]
        genome.models[[kk]][[k]] <- as.matrix(genome.models1[ , ptr.models[[k]] ])
        betatXtY.2[[kk]][[k]] <- as.matrix(betatXtY.22[ ptr.models[[k]], 2, drop=FALSE ])
        betatXtY.3[[kk]][[k]] <- as.matrix(betatXtY.33[ ptr.models[[k]], 2, drop=FALSE ])
        model.names[[kk]][[k]] <- betatXtY.22[ ptr.models[[k]], 1 ]

        R2.indiv2 <- 2*log(abs( betatXtY.2[[kk]][[k]] )) - log(diag(Sigma.prs[[kk]][[k]]))
        R2.indiv3 <- 2*log(abs( betatXtY.3[[kk]][[k]] )) - log(diag(Sigma.prs[[kk]][[k]]))
        R2.2 <- exp(R2.indiv2 - max(R2.indiv2,na.rm=TRUE))
        R2.3 <- exp(R2.indiv3 - max(R2.indiv3,na.rm=TRUE))

        ptr.use[[k]] <- which( diag(Sigma.prs[[kk]][[k]])!=0 &
                               R2.2 / max(R2.2[betatXtY.2[[kk]][[k]]>0],na.rm=TRUE) > 0.2 &
                               R2.3 / max(R2.3[betatXtY.3[[kk]][[k]]>0],na.rm=TRUE) > 0.2 &
                               betatXtY.2[[kk]][[k]] > 0 & betatXtY.3[[kk]][[k]] > 0 )
        n.prs <- length(ptr.use[[k]])

        if( n.prs>0 ){
            Sigma.prs[[kk]][[k]] <- Sigma.prs[[kk]][[k]][ ptr.use[[k]], ptr.use[[k]], drop=FALSE ]
            genome.models[[kk]][[k]] <- genome.models[[kk]][[k]][,ptr.use[[k]], drop=FALSE ]
            betatXtY.2[[kk]][[k]] <- betatXtY.2[[kk]][[k]][ptr.use[[k]],,drop=FALSE]
            betatXtY.3[[kk]][[k]] <- betatXtY.3[[kk]][[k]][ptr.use[[k]],,drop=FALSE]
            model.names[[kk]][[k]] <- model.names[[kk]][[k]][ptr.use[[k]],1]

            prs.norm <- 1 / sqrt(diag(Sigma.prs[[kk]][[k]]))
            genome.models[[kk]][[k]] <- genome.models[[kk]][[k]] %*% diag(prs.norm,nrow=n.prs)
            Sigma.prs[[kk]][[k]] <- diag(prs.norm,nrow=n.prs) %*% Sigma.prs[[kk]][[k]] %*% diag(prs.norm,nrow=n.prs)
            betatXtY.2[[kk]][[k]] <- betatXtY.2[[kk]][[k]] * prs.norm
            betatXtY.3[[kk]][[k]] <- betatXtY.3[[kk]][[k]] * prs.norm

            lambda.range <- c( 1, 10*n.test )
            lambda.opt <- optimise( enet.R2, interval=lambda.range,
                                   n=n.test, Sigma.prs=Sigma.prs[[kk]][[k]],
                                   betatXtY.2=betatXtY.2[[kk]][[k]], betatXtY.3=betatXtY.3[[kk]][[k]],
                                   maximum=TRUE )
            lambda[kk,k] <- lambda.opt$maximum

            prs.weights <- solve( diag(lambda[kk,k],n.prs) + n.test*Sigma.prs[[kk]][[k]] ) %*% betatXtY.2[[kk]][[k]]
            R2.model[kk,k] <- 2*log(abs(t(prs.weights) %*% betatXtY.3[[kk]][[k]])) -
                log( diag(t(prs.weights) %*% Sigma.prs[[kk]][[k]] %*% prs.weights) )
        }
    }
    R2.indiv[[kk]] <- 2*log(abs( betatXtY.2[[kk]][[n.models]] + betatXtY.3[[kk]][[n.models]] ))
    all.model.names <- union( all.model.names, model.names[[kk]][[n.models]]$V1 )
}

R2.indiv1 <- rep( 0, length(all.model.names) )
count.indiv1 <- rep( 0, length(all.model.names) )
names(R2.indiv1) <- all.model.names
for( kk in 1:opt$n.folds ){
    ptr <- match( model.names[[kk]][[n.models]]$V1, all.model.names )
    R2.indiv1[ptr] <- R2.indiv1[ptr] + R2.indiv[[kk]]
    count.indiv1[ptr] <- count.indiv1[ptr] + 1
}
s2 <- order( R2.indiv1, decreasing=TRUE )
best.model.name <- all.model.names[s2[1]]

k <- get_mode( apply( R2.model, 1, order, decreasing=TRUE )[1,] )
lambda1 <- median( lambda[,k], na.rm=TRUE )
print( c( k, lambda1 ) )
ensembl.model <- as.data.frame(matrix( ncol=4, nrow=length(snps), data=0 ))
best.model <- as.data.frame(matrix( ncol=4, nrow=length(snps), data=0 ))
colnames(ensembl.model) <- c('SNP.id', 'ALLELE1', 'ALLELE0', 'beta' )
colnames(best.model) <- c('SNP.id', 'ALLELE1', 'ALLELE0', paste0("beta_",best.model.name) )
for( kk in 1:opt$n.folds ){
    ptr <- match( alleles[[kk]]$genome.models1.rownames, snps )
    ensembl.model[ptr,1:3] <- alleles[[kk]][,1:3]
    best.model[ptr,1:3] <- alleles[[kk]][,1:3]

    ptr.model <- match( best.model.name, model.names[[kk]][[n.models]]$V1 )
    if( !is.na(ptr.model) ){
        best.model[ptr,4] <- best.model[ptr,4] + genome.models[[kk]][[n.models]][,ptr.model]
    }
    n.prs <- nrow(Sigma.prs[[kk]][[k]])
    prs.weights <- solve( diag(lambda1,n.prs) + n.test.all*Sigma.prs[[kk]][[k]] ) %*%
        ( betatXtY.2[[kk]][[k]] + betatXtY.3[[kk]][[k]] )
    ensembl.model[ptr,4] <- ensembl.model[ptr,4] + genome.models[[kk]][[k]] %*% as.matrix( prs.weights )
}
fwrite( ensembl.model,
       paste0(opt$workdir,'/snp_weights_weighted_model.dat'),
       row.names=FALSE, col.names=TRUE, sep=" " )
fwrite( best.model,
       paste0(opt$workdir,'/snp_weights_best_model.dat'),
       row.names=FALSE, col.names=TRUE, sep=" " )


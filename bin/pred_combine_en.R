library(data.table)
library(boot)
library(glmnet)
library(doMC)
library("optparse")
source('~/BridgePRS/bin/functions.R')

var.explained <- function(data,ptr){
    ptr.X <- grep('X',colnames(data))
    ptr.PRS <- grep('PRS',colnames(data))
    fit0 <- summary(lm( data$y ~ 0 + as.matrix(data[,ptr.X]), subset=ptr ))
    fit1 <- summary(lm( data$y ~ 0 + as.matrix(data[,ptr.X]) + as.matrix(data[,ptr.PRS]), subset=ptr ))
    R2 <- 1 - (1-fit1$adj.r.squared) / (1-fit0$adj.r.squared)
    return(R2)
}

option_list = list(
    make_option(c("--pred.dir1"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--pred.dir2"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--pred1"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--pred2"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--outdir"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--pop2"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--valid.data"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--test.data"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--cov.names"), type="character", default="000",
                help="", metavar="character"),
    make_option(c("--pheno.name"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--ids.col"), type="logical", default=TRUE,
                help="", metavar="logical"),
    make_option(c("--n.cores"), type="numeric", default=1,
                help="Number of processors for mclapply to use", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print(opt)

tmp <- t(data.frame(opt))
rownames(tmp) <- names(opt)
logfile <- paste0(opt$outdir,"/",opt$pop2,"_weighted_combined_var_explained.log")
write.table(tmp,file=logfile,quote=FALSE,col.names=FALSE)

if( opt$cov.names!="000" ){
    cov.names <- strsplit( opt$cov.names, ',' )[[1]]
}
registerDoMC(cores = opt$n.cores)

#### Prediction ####

pred.dir1 <- ifelse( is.null(opt$pred.dir1), opt$outdir, opt$pred.dir1 )
pred.dir2 <- ifelse( is.null(opt$pred.dir2), opt$outdir, opt$pred.dir2 )

pred1 <- fread(paste0(pred.dir1, "/", opt$pop,'_',opt$pred1,'_all_preds_test.dat'),
               data.table=FALSE)
target <- fread(opt$test.data,data.table=FALSE)
if( opt$ids.col ){
    target <- target[match(pred1$id,target$IID),]
    pred1 <- pred1[,-1]
}
v.target <- var(target[,opt$pheno.name])

nfolds <- ifelse( nrow(target)<2000, nrow(target), 50 )

write( paste("Test data of", nrow(target), "samples"), file=logfile )

fit.ridge1 <- cv.glmnet( y=target[,opt$pheno.name], x=as.matrix(pred1),
                        family="gaussian",
                        alpha=0, parallel=TRUE, nfolds=nfolds, grouped=FALSE )

if( !is.null(opt$pred2) ){
    pred2 <- fread(paste0(pred.dir2, "/", opt$pop,'_',opt$pred2,'_all_preds_test.dat'),
                   data.table=FALSE)
    if( opt$ids.col ){
        pred2 <- pred2[,-1]
    }
    X <- as.matrix(cbind( pred1, pred2 ))
    fit.ridge <- cv.glmnet( y=target[,opt$pheno.name], x=X, family="gaussian",
                           alpha=0, parallel=TRUE, nfolds=nfolds, grouped=FALSE )
    fit.ridge2 <- cv.glmnet( y=target[,opt$pheno.name], x=as.matrix(pred2),
                            family="gaussian",
                            alpha=0, parallel=TRUE, nfolds=nfolds, grouped=FALSE )

    w.ridge <- getGlmnetFit( fit.ridge, X, s='lambda.min', sparse=FALSE )[-1]
    w.ridge1 <- getGlmnetFit( fit.ridge1, as.matrix(pred1), s='lambda.min', sparse=FALSE )[-1]
    w.ridge2 <- getGlmnetFit( fit.ridge2, as.matrix(pred2), s='lambda.min', sparse=FALSE )[-1]

    n <- fit.ridge$glmnet.fit$nobs
    ptr.min <- which(fit.ridge$lambda==fit.ridge$lambda.min)
    ptr.min1 <- which(fit.ridge1$lambda==fit.ridge1$lambda.min)
    ptr.min2 <- which(fit.ridge2$lambda==fit.ridge2$lambda.min)
    mse <- c( fit.ridge$cvm[ptr.min], fit.ridge1$cvm[ptr.min1], fit.ridge2$cvm[ptr.min2] )
#logL2 <- -n*( mse - as.numeric(mse[1]) ) / (2*v.target)
#probM2 <- exp(logL2) / sum(exp(logL2))
    logL <- -n*log(mse) / 2
    logL <- logL - min(logL)
    probM.ridge <- exp(logL) / sum(exp(logL))
}


#### Validation ####
if( opt$valid.data!=0 ){
    pred1 <- fread(paste0(pred.dir1, "/", opt$pop,'_',opt$pred1,
                          '_all_preds_valid.dat'), data.table=FALSE)

    target <- fread(opt$valid.data,data.table=FALSE)
    if( opt$ids.col ){
        target <- target[match(pred1$id,target$IID),]
        pred1 <- pred1[,-1]
    }
    write( paste("Validation data of", nrow(target), "samples"), file=logfile )

    if( opt$cov.names!="000" ){
        covs <-  cbind( 1, target[,cov.names] )
    }else{
        covs <- matrix(ncol=1,nrow=nrow(target),data=1)
    }
    colnames(covs) <- paste('X',colnames(covs),sep='.')
    prs.ridge1 <- predict( fit.ridge1, as.matrix(pred1), s='lambda.min' )

    data <- data.frame( target[,opt$pheno.name], prs.ridge1, covs )
    colnames(data)[1:2] <- c('y','PRS')
    b <- boot(data,var.explained,stype="i",R=10000,parallel='multicore',ncpus=opt$n.cores)
    ci <- boot.ci(b,type='norm')
    ptr.min <- which(fit.ridge1$lambda==fit.ridge1$lambda.min)
    VE.ridge1 <- c( b$t0,ci$normal[-1], fit.ridge1$cvm[ptr.min], fit.ridge1$cvsd[ptr.min] )
    pred1 <- fread(paste0(pred.dir1, "/", opt$pop,'_',opt$pred1,
                          '_all_preds_valid.dat'), data.table=FALSE)

    ptr <- c( 1, grep("_10_",colnames(pred1)), grep("_15_",colnames(pred1)),
             grep("_20_",colnames(pred1)), grep("_50_",colnames(pred1)) )
    pred1 <- pred1[,ptr]

    target <- fread(opt$valid.data,data.table=FALSE)
    if( opt$ids.col ){
        target <- target[match(pred1$id,target$IID),]
        pred1 <- pred1[,-1]
    }
    write( paste("Validation data of", nrow(target), "samples"), file=logfile )

    if( opt$cov.names!="000" ){
        covs <-  cbind( 1, target[,cov.names] )
    }else{
        covs <- matrix(ncol=1,nrow=nrow(target),data=1)
    }
    colnames(covs) <- paste('X',colnames(covs),sep='.')
    prs.ridge1 <- predict( fit.ridge1, as.matrix(pred1), s='lambda.min' )

    data <- data.frame( target[,opt$pheno.name], prs.ridge1, covs )
    colnames(data)[1:2] <- c('y','PRS')
    b <- boot(data,var.explained,stype="i",R=10000,parallel='multicore',ncpus=opt$n.cores)
    ci <- boot.ci(b,type='norm')
    ptr.min <- which(fit.ridge1$lambda==fit.ridge1$lambda.min)
    VE.ridge1 <- c( b$t0,ci$normal[-1], fit.ridge1$cvm[ptr.min], fit.ridge1$cvsd[ptr.min] )

    if( !is.null(opt$pred2) ){
        pred2 <- fread(paste0( pred.dir2, "/", opt$pop,'_',opt$pred2,
                              '_all_preds_valid.dat'), data.table=FALSE)
        if( opt$ids.col ){
            pred2 <- pred2[,-1]
        }
        X <- as.matrix(cbind( pred1, pred2 ) )
        prs.ridge <- predict( fit.ridge, X, s='lambda.min' )
        prs.ridge2 <- predict( fit.ridge2, as.matrix(pred2), s='lambda.min' )

        ptr.fit1 <- match( names(w.ridge1), names(w.ridge) )
        ptr.fit2 <- match( names(w.ridge2), names(w.ridge) )

        data <- data.frame( target[,opt$pheno.name], prs.ridge, covs )
        colnames(data)[1:2] <- c('y','PRS')
        b <- boot(data,var.explained,stype="i",R=10000,parallel='multicore',ncpus=opt$n.cores)
        ci <- boot.ci(b,type='norm')
        ptr.min <- which(fit.ridge$lambda==fit.ridge$lambda.min)
        VE.ridge <- c(b$t0,ci$normal[-1], fit.ridge$cvm[ptr.min], fit.ridge$cvsd[ptr.min] )

        data <- data.frame( target[,opt$pheno.name], prs.ridge2, covs )
        colnames(data)[1:2] <- c('y','PRS')
        b <- boot(data,var.explained,stype="i",R=10000,parallel='multicore',ncpus=opt$n.cores)
        ci <- boot.ci(b,type='norm')
        ptr.min <- which(fit.ridge2$lambda==fit.ridge2$lambda.min)
        VE.ridge2 <- c(b$t0,ci$normal[-1], fit.ridge2$cvm[ptr.min], fit.ridge2$cvsd[ptr.min] )

        prs.weighted <- apply( cbind( prs.ridge, prs.ridge1, prs.ridge2 ) %*% diag(probM.ridge), 1, sum )
        data <- data.frame( target[,opt$pheno.name], prs.weighted, covs )
        colnames(data)[1:2] <- c('y','PRS')
        b <- boot(data,var.explained,stype="i",R=10000,parallel='multicore',ncpus=opt$n.cores)
        ci <- boot.ci(b,type='norm')
        ptr.min <- which(fit.ridge1$lambda==fit.ridge1$lambda.min)
        VE.ridge.w <- c(b$t0,ci$normal[-1], fit.ridge1$cvm[ptr.min], fit.ridge1$cvsd[ptr.min] )
        write.table( data.frame( target$IID, prs.weighted ),
                    paste0(opt$outdir,"/",opt$pop2,"_weighted_combined_preds.dat"),
                    col.names=FALSE, row.names=FALSE, quote=FALSE )

        out <- rbind( VE.ridge, VE.ridge1, VE.ridge2, VE.ridge.w )
        out <- cbind( c(probM.ridge,0), out )
        colnames(out) <- c('Prob','Est','2.5%','97.5%','cv.dev','cv.dev.sd')
        rownames(out) <- c('Ridge','Ridge1','Ridge2','Ridge.w')
        print(out)
        write.csv( out, paste0(opt$outdir,"/",opt$pop2,"_weighted_combined_var_explained.txt"), row.names=TRUE )
    }
    if( is.null(opt$pred2) ){
        out <- rbind( VE.ridge1 )
        colnames(out) <- c('Est','2.5%','97.5%','cv.dev','cv.dev.sd')
        rownames(out) <- c('Ridge')
        print(out)
        write.csv( out, paste0(opt$outdir,"/",opt$pop2,"_weighted_combined_var_explained.txt"), row.names=TRUE )
    }
}

tmp <- strsplit( names(w.ridge1), "_" )
tau <- as.numeric(sapply(tmp,getElement,2))
tau.weights2 <- tapply( w.ridge1, tau, sum )

tmp <- strsplit( names(w.ridge2), "_" )
lambda <- as.numeric(sapply(tmp,getElement,2))
alpha <- as.numeric(sapply(tmp,getElement,3))
lambda.weights2 <- tapply( w.ridge2, lambda, sum )
alpha.weights2 <- tapply( w.ridge2, alpha, sum )

tmp <- strsplit( names(w.ridge), "_" )
l <- sapply(tmp,length)
lambda <- as.numeric(sapply(tmp,getElement,2))
alpha <- as.numeric(sapply(tmp,getElement,3))
tau <- as.numeric(sapply(tmp,getElement,2))
tau.weights1 <- tapply( w.ridge[l==3], tau[l==3], sum )
alpha.weights1 <- tapply( w.ridge[l==4], alpha[l==4], sum )
lambda.weights1 <- tapply( w.ridge[l==4], lambda[l==4], sum )

tau.weights <- tau.weights2 * probM.ridge[2] + tau.weights1 * probM.ridge[1]
alpha.weights <- alpha.weights2 * probM.ridge[3] + alpha.weights1 * probM.ridge[1]
lambda.weights <- lambda.weights2 * probM.ridge[3] + lambda.weights1 * probM.ridge[1]

tau.weights <- cbind( tau.weights1, tau.weights2, tau.weights )
alpha.weights <- cbind( alpha.weights1, alpha.weights2, alpha.weights )
lambda.weights <- cbind( lambda.weights1, lambda.weights2, lambda.weights )

write.table( tau.weights, paste0(opt$outdir,"/",opt$pop2,"_tau_weights.dat"), row.names=TRUE )
write.table( alpha.weights, paste0(opt$outdir,"/",opt$pop2,"_alpha_weights.dat"), row.names=TRUE )
write.table( lambda.weights, paste0(opt$outdir,"/",opt$pop2,"_lambda_weights.dat"), row.names=TRUE )

# Contribution of stage 2 model
w.ridge11 <- w.ridge[1:length(w.ridge1)] * probM.ridge[1] +
    w.ridge1 * probM.ridge[2]
# Contribution of single ancestry model
w.ridge22 <- w.ridge[-(1:length(w.ridge1))] * probM.ridge[1] +
    w.ridge2 * probM.ridge[3]

beta.bar.genome <- as.data.frame(matrix( nrow=0, ncol=4 ))
for( chr in 1:22 ){
    beta.bar1 <- fread(paste0(pred.dir1,"/models/",opt$pop,"_",
                              opt$pred1,"_beta_bar_chr",chr,".txt.gz"),
                       data.table=FALSE)
    beta.bar2 <- fread(paste0(pred.dir2,"/models/",opt$pop,"_",
                              opt$pred2,"_beta_bar_chr",chr,".txt.gz"),
                       data.table=FALSE)
    beta.bar1 <- beta.list(beta.bar1)
    beta.bar2 <- beta.list(beta.bar2)

    kl <- fread(paste0(pred.dir1,"/models/",opt$pop,"_",
                       opt$pred1,"_KLdist_chr",chr,".txt.gz"),
                data.table=FALSE)[,-1]
    if( chr==1 ){
        p.clump <- sapply(sapply(beta.bar1,getElement,'p.value'),getElement,1)
        F <- ecdf(p.clump)
        q <- F(10^(-(1:8)))
        kl.thresh <- apply( kl, 2, quantile, 1-q, na.rm=TRUE )
    }

    tmp <- strsplit( names(w.ridge11), '_' )
    w.thresh <- sapply(tmp,getElement,3)
    w.beta.bar <- paste( sapply(tmp,getElement,1),
                        sapply(tmp,getElement,2), sep='_' )
    u.w.thresh <- unique(w.thresh)
    beta.bar11 <- as.data.frame(matrix( nrow=0, ncol=4 ))
    for( i in 1:length(beta.bar1) ){
        beta.bar <- as.data.frame(matrix( nrow=nrow(beta.bar1[[i]]), ncol=4 ))
        beta.bar[,1:3] <- beta.bar1[[i]][,c('snp','effect.allele','ref.allele')]
        beta.bar[,4] <- 0
        for( j in 1:length(u.w.thresh) ){
            ptr.weight <- which( w.thresh==u.w.thresh[j] )
            ptr.beta.bar <- match( w.beta.bar[ptr.weight], colnames(beta.bar1[[i]]) )
            ptr.kl <- match( colnames(beta.bar1[[i]])[ptr.beta.bar], colnames(kl) )
            ptr <- which( kl[i,ptr.kl] > kl.thresh[j,ptr.kl] )
            tmp <- matrix( nrow=nrow(beta.bar1[[i]]), ncol=length(ptr.beta.bar), data=0 )
            tmp[,ptr] <- as.matrix(beta.bar1[[i]][,ptr.beta.bar[ptr]])
            beta.bar[,4] <- beta.bar[,4] + tmp %*% w.ridge11[ptr.weight]
        }
        if( sum(beta.bar[,4])!=0 ){
            beta.bar11 <- rbind( beta.bar11, beta.bar )
        }
    }

    tmp <- strsplit( names(w.ridge22), '_' )
    w.thresh <- as.numeric(sapply(tmp,getElement,4))
    w.beta.bar <- paste( sapply(tmp,getElement,1),
                        sapply(tmp,getElement,2),
                        sapply(tmp,getElement,3), sep='_' )
    u.w.thresh <- unique(w.thresh)
    p.clump <- sapply(sapply(beta.bar2,getElement,'p.value'),getElement,1)
    beta.bar22 <- as.data.frame(matrix( nrow=0, ncol=4 ))
    for( i in 1:length(beta.bar2) ){
        j <- 1
        beta.bar <- as.data.frame(matrix( nrow=nrow(beta.bar2[[i]]), ncol=4 ))
        beta.bar[,1:3] <- beta.bar2[[i]][,c('snp','effect.allele','ref.allele')]
        beta.bar[,4] <- 0
        while( j<=length(u.w.thresh) & p.clump[i] < u.w.thresh[j] ){
            ptr.weight <- which( w.thresh==u.w.thresh[j] )
            ptr.beta.bar <- match( w.beta.bar[ptr.weight], names(beta.bar2[[i]]) )
            beta.bar[,4] <- beta.bar[,4] +
                as.matrix(beta.bar2[[i]][,ptr.beta.bar]) %*% w.ridge22[ptr.weight]
            j <- j+1
        }
        if( j>1 ){
            beta.bar22 <- rbind( beta.bar22, beta.bar )
        }
    }
    i.snps <- intersect(beta.bar11$V1,beta.bar22$V1)
    ptr1 <- match( i.snps, beta.bar11$V1 )
    ptr2 <- match( i.snps, beta.bar22$V1 )
    beta.bar11[ ptr1, 4] <- beta.bar11[ ptr1, 4] + beta.bar22[ ptr2, 4]
    beta.bar.chr <- rbind( beta.bar11, beta.bar22[-ptr2,] )
    beta.bar.genome <- rbind( beta.bar.genome, beta.bar.chr )
}
write.table( beta.bar.genome,
          paste0(opt$outdir,"/",opt$pop2,"_weighted_combined_snp_weights.dat"),
          col.names=TRUE, row.names=FALSE, quote=FALSE )

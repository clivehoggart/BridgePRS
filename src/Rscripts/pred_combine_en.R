library(data.table)
library(boot)
library(glmnet)
library(doMC)
library("optparse")

calc_ve <- function(prs, myfit) {
    data <- data.frame(target[, opt$pheno.name], prs, covs)
    colnames(data)[1:2] <- c("y", "PRS")
    b <- boot(data, function(x, ptr) var_explained(x, ptr, family),
              stype = "i", R = 10000, parallel = "multicore",
              ncpus = opt$n.cores)
    ci <- boot.ci(b, type = "norm")

    cvm <- 1
    cvsd <- 1
    if( !is.null(myfit) ){
        ptr.min <- which(myfit$lambda == myfit$lambda.min)
        cvm <-  myfit$cvm[ptr.min]
        cvsd <- myfit$cvsd[ptr.min]
    }

    c(b$t0, ci$normal[-1], cvm, cvsd)
}

var_explained <- function(data, ptr, fam = "gaussian") {
  # Covariates
  ptr_x <- grep("X", colnames(data))
  # PRS
  ptr_prs <- grep("PRS", colnames(data))

  # Difference in R2
  if (fam == "gaussian") {
    # Null Model
    fit0 <- lm(data$y ~ 0 + as.matrix(data[, ptr_x]), subset = ptr)
    # Full Model
    fit1 <- lm(
      data$y ~ 0 + as.matrix(data[, ptr_x]) + as.matrix(data[, ptr_prs]),
      subset = ptr)
    r2 <- (1 - (1 - summary(fit1)$adj.r.squared) /
               (1 - summary(fit0)$adj.r.squared))
  } else {
    if (fam == "binomial" && all(range(data$y) == c(1, 2))) {
      data$y <- data$y - 1
    }
    # Null Model
    fit0 <- glm(data$y ~ 0 + as.matrix(data[, ptr_x]),
                subset = ptr, family = fam)
    # Full Model
    fit1 <- glm(
      data$y ~ 0 + as.matrix(data[, ptr_x]) + as.matrix(data[, ptr_prs]),
      subset = ptr, family = fam)
    n <- length(fit0$y)
    d_null <- -2 * logLik(fit0)[1]
    d_full <- -2 * logLik(fit1)[1]
    r2 <- (1 - exp((d_full - d_null) / n)) / (1 - exp(-d_null / n))
  }
  return(r2)
}

option_list = list(
    make_option(c("--fpath"), type="character", default=NULL,
                help="Function File Path", metavar="character"),
    make_option(c("--pred1"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--pred2"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--models1"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--models2"), type="character", default=NULL,
                help="", metavar="character"),
    make_option(c("--outfile"), type="character", default=NULL,
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
                help="Number of processors for mclapply to use", metavar="character"),
    make_option(c("--binary"), type="numeric", default=0,
                help="Indicator for binary outcome", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print(opt)
source(opt$fpath)

if( opt$binary==1 ){
    family <- "binomial"
}else{
    family <- "gaussian"
}

tmp <- t(data.frame(opt))
rownames(tmp) <- names(opt)


# TADE EDIT - Renaming LogFile For Pred1/2

#logfile <- paste0(opt$outfile,"_weighted_combined_var_explained.log")
if( is.null(opt$pred2) ){
    logfile <- paste0(opt$outfile,"_var_explained.log")
    wtfile  <- paste0(opt$outfile,"_snp_weights.dat")

}else{
    logfile <- paste0(opt$outfile,"_weighted_combined_var_explained.log")
    wtfile <- paste0(opt$outfile,"_weighted_combined_snp_weights.dat")
}





write.table(tmp,file=logfile,quote=FALSE,col.names=FALSE)
if( opt$cov.names!="000" ){
    cov.names <- strsplit( opt$cov.names, ',' )[[1]]
}
registerDoMC(cores = opt$n.cores)

#### Prediction ####

pred1 <- fread(paste0(opt$pred1,'_all_preds_test.dat'), data.table=FALSE)
target <- fread(opt$test.data,data.table=FALSE)
if( opt$ids.col ){
    target <- target[match(pred1$id,target$IID),]
    pred1 <- pred1[,-1]
}
v.target <- var(target[,opt$pheno.name])

#ptr <- c( grep("_0.5_", colnames(pred1)),
#         grep("_1_", colnames(pred1)),
#         grep("_2_", colnames(pred1)),
#         grep("_5_", colnames(pred1)),
#         grep("_10_", colnames(pred1)) )
#pred1 <- pred1[,ptr]

nfolds <- ifelse( nrow(target)<500, nrow(target), 500 )

write( paste("Test data of", nrow(target), "samples"), file=logfile )

fit.ridge1 <- cv.glmnet( y=target[,opt$pheno.name], x=as.matrix(pred1),
                        family=family,
                        alpha=0, parallel=TRUE, nfolds=nfolds, grouped=FALSE )
#w.ridge1 <- getGlmnetFit( fit.ridge1, as.matrix(pred1), s='lambda.min', sparse=FALSE )[-1]
w.ridge1 <- as.vector(coef( fit.ridge1, s='lambda.min' ))[-1]
names(w.ridge1) <- colnames(pred1)

if( !is.null(opt$pred2) ){
    pred2 <- fread(paste0(opt$pred2,'_all_preds_test.dat'), data.table=FALSE)
    if( opt$ids.col ){
        pred2 <- pred2[,-1]
    }
    X <- as.matrix(cbind( pred1, pred2 ))
    fit.ridge <- cv.glmnet( y=target[,opt$pheno.name], x=X, family=family,
                           alpha=0, parallel=TRUE, nfolds=nfolds, grouped=FALSE )
    fit.ridge2 <- cv.glmnet( y=target[,opt$pheno.name], x=as.matrix(pred2),
                            family=family,
                            alpha=0, parallel=TRUE, nfolds=nfolds, grouped=FALSE )

#    w.ridge <- getGlmnetFit( fit.ridge, X, s='lambda.min', sparse=FALSE )[-1]
#    w.ridge2 <- getGlmnetFit( fit.ridge2, as.matrix(pred2), s='lambda.min', sparse=FALSE )[-1]
    w.ridge <- as.vector(coef( fit.ridge, s='lambda.min' )[-1])
    w.ridge2 <- as.vector(coef( fit.ridge2, s='lambda.min' )[-1])
    names(w.ridge) <- colnames(X)
    names(w.ridge2) <- colnames(pred2)

    n <- fit.ridge$glmnet.fit$nobs
    ptr.min <- which(fit.ridge$lambda==fit.ridge$lambda.min)
    ptr.min1 <- which(fit.ridge1$lambda==fit.ridge1$lambda.min)
    ptr.min2 <- which(fit.ridge2$lambda==fit.ridge2$lambda.min)
    mse <- c( fit.ridge$cvm[ptr.min], fit.ridge1$cvm[ptr.min1], fit.ridge2$cvm[ptr.min2] )
#logL2 <- -n*( mse - as.numeric(mse[1]) ) / (2*v.target)
#probM2 <- exp(logL2) / sum(exp(logL2))
    logL <- -n*log(mse) / 2
    logL <- logL - max(logL)
    probM.ridge <- exp(logL) / sum(exp(logL))
}

# alpha

#### Validation ####
if( opt$valid.data!=0 ){
    pred1 <- fread(paste0(opt$pred1,'_all_preds_valid.dat'), data.table=FALSE)
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
    VE.ridge1 <- calc_ve(prs.ridge1, fit.ridge1)

    if( !is.null(opt$pred2) ){
        pred2 <- fread(paste0( opt$pred2, '_all_preds_valid.dat'), data.table=FALSE)
        if( opt$ids.col ){
            pred2 <- pred2[,-1]
        }
        X <- as.matrix(cbind( pred1, pred2 ) )
        prs.ridge <- predict( fit.ridge, X, s='lambda.min' )
        prs.ridge2 <- predict( fit.ridge2, as.matrix(pred2), s='lambda.min' )

        ptr.fit1 <- match( names(w.ridge1), names(w.ridge) )
        ptr.fit2 <- match( names(w.ridge2), names(w.ridge) )

        VE.ridge <- calc_ve(prs.ridge, fit.ridge)
        VE.ridge2 <- calc_ve(prs.ridge2, fit.ridge2)

        prs.weighted <- apply( cbind( prs.ridge, prs.ridge1, prs.ridge2 ) %*% diag(probM.ridge), 1, sum )
        VE.ridge.w <- calc_ve(prs.weighted, NULL)

        # TADE - NEED THE PHENOTYPE AND COLUMN NAMES FOR DOWNSTREAM ANALYSIS #
        out <- data.frame( target$IID, target[,opt$pheno.name], prs.ridge, prs.ridge1, prs.ridge2, prs.weighted )
        colnames(out) <- c('---','pheno','prs','ridge1','ridge2','weighted')

        # TADE - NEED THE COLUMN NAMES FOR DOWNSTREAM
        #write.table( out,paste0(opt$outfile,"_weighted_combined_preds.dat"),col.names=FALSE, row.names=FALSE, quote=FALSE )

        #my_data <- data.frame( target$IID, target[,opt$pheno.name], prs.ridge1)
        write.table( out,paste0(opt$outfile,"_weighted_combined_preds.dat"),row.names=FALSE, quote=FALSE )

        out <- rbind( VE.ridge, VE.ridge1, VE.ridge2, VE.ridge.w )
        out <- cbind( c(probM.ridge,1), out )
        colnames(out) <- c('Prob','Est','2.5%','97.5%','cv.dev','cv.dev.sd')
        rownames(out) <- c('Ridge','Ridge1','Ridge2','Ridge.w')
        print(out)
        write.csv( out, paste0(opt$outfile,"_weighted_combined_var_explained.txt"), row.names=TRUE )
    }
    if( is.null(opt$pred2) ){
        out <- rbind( VE.ridge1 )
        colnames(out) <- c('Est','2.5%','97.5%','cv.dev','cv.dev.sd')
        rownames(out) <- c('Ridge')
        print(out)

        # TADE - Renaming Again (not weighted/combined)
        write.csv( out, paste0(opt$outfile,"_var_explained.txt"), row.names=TRUE )

        # TADE - Adding Single Stage Preds for single stage #
        data <- data.frame( target[  ,opt$pheno.name], prs.ridge1, covs )
        my_data <- data.frame( target$IID, target[,opt$pheno.name], prs.ridge1)
        colnames(my_data) <- c('---','pheno','prs')
        write.table( my_data, paste0(opt$outfile,"_preds.dat"), row.names = FALSE, quote=FALSE)

    }
}

tmp <- strsplit( names(w.ridge1), "_" )
lambda <- as.numeric(sapply(tmp,getElement,2))
alpha <- as.numeric(sapply(tmp,getElement,3))
lambda.weights2 <- tapply( w.ridge1, lambda, sum )
alpha.weights2 <- tapply( w.ridge1, alpha, sum )
if( is.null(opt$pred2) ){
    #write.table( alpha.weights2,  paste0(opt$outdir,"_alpha_weights.dat"), row.names=TRUE )
    #write.table( lambda.weights2, paste0(opt$outdir,"_lambda_weights.dat"), row.names=TRUE )
    write.table( alpha.weights2,  paste0(opt$outfile,"_alpha_weights.dat"), row.names=TRUE )
    write.table( lambda.weights2, paste0(opt$outfile,"_lambda_weights.dat"), row.names=TRUE )
}

if( !is.null(opt$pred2) ){
    tmp <- strsplit( names(w.ridge2), "_" )
    tau <- as.numeric(sapply(tmp,getElement,2))
    tau.weights2 <- tapply( w.ridge2, tau, sum )

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

    write.table( tau.weights, paste0(opt$outfile,"_tau_weights.dat"), row.names=TRUE )
    write.table( alpha.weights, paste0(opt$outfile,"_alpha_weights.dat"), row.names=TRUE )
    write.table( lambda.weights, paste0(opt$outfile,"_lambda_weights.dat"), row.names=TRUE )
}

w.ridge11 <- w.ridge1
if( !is.null(opt$pred2) ){
# Contribution of single ancestry model
    w.ridge11 <- w.ridge[1:length(w.ridge1)] * probM.ridge[1] + w.ridge1 * probM.ridge[2]
# Contribution of stage 2 model
    w.ridge22 <- w.ridge[-(1:length(w.ridge1))] * probM.ridge[1] + w.ridge2 * probM.ridge[3]
}


#print('hello')
#match
beta.bar.genome <- as.data.frame(matrix( nrow=0, ncol=4 ))
for( chr in 1:22 ){
    beta.bar1 <- fread(paste0(opt$models1,"_beta_bar_chr",chr,".txt.gz"),
                       data.table=FALSE)
    beta.bar1 <- beta.list(beta.bar1)

    if( !is.null(opt$models2) ){
        beta.bar2 <- fread(paste0(opt$models2,"_beta_bar_chr",chr,".txt.gz"),data.table=FALSE)
        beta.bar2 <- beta.list(beta.bar2)

        kl <- fread(paste0(opt$models2,"_KLdist_chr",chr,".txt.gz"),
                    data.table=FALSE)[,-1]
        if( chr==1 ){
            p.clump <- sapply(sapply(beta.bar2,getElement,'p.value'),getElement,1)
            F <- ecdf(p.clump)
            q <- F(10^(-(1:8)))
            kl.thresh <- apply( kl, 2, quantile, 1-q, na.rm=TRUE )
        }

        tmp <- strsplit( names(w.ridge22), '_' )
        w.thresh <- sapply(tmp,getElement,3)
        w.beta.bar <- paste( sapply(tmp,getElement,1),sapply(tmp,getElement,2), sep='_' )
        u.w.thresh <- unique(w.thresh)
        beta.bar22 <- as.data.frame(matrix( nrow=0, ncol=4 ))

        for( i in 1:length(beta.bar2) ){
            beta.bar <- as.data.frame(matrix( nrow=nrow(beta.bar2[[i]]), ncol=4 ))
            beta.bar[,1:3] <- beta.bar2[[i]][,c('snp','effect.allele','ref.allele')]
            beta.bar[,4] <- 0
            for( j in 1:length(u.w.thresh) ){
                ptr.weight <- which( w.thresh==u.w.thresh[j] )
                ptr.beta.bar <- match( w.beta.bar[ptr.weight], colnames(beta.bar2[[i]]) )
                ptr.kl <- match( colnames(beta.bar2[[i]])[ptr.beta.bar], colnames(kl) )
                ptr <- which( kl[i,ptr.kl] > kl.thresh[j,ptr.kl] )
                tmp <- matrix( nrow=nrow(beta.bar2[[i]]), ncol=length(ptr.beta.bar), data=0 )
                tmp[,ptr] <- as.matrix(beta.bar2[[i]][,ptr.beta.bar[ptr]])
                beta.bar[,4] <- beta.bar[,4] + tmp %*% w.ridge22[ptr.weight]
            }
            if( sum(beta.bar[,4])!=0 ){
                beta.bar22 <- rbind( beta.bar22, beta.bar )
            }
        }
    }

    tmp <- strsplit( names(w.ridge11), '_' )

    ## NOTE - THIS CAUSES AN ERROR WHEN RAN WITH SINGLE PRIOR RESULT - tmp, getElement4 ##
    w.thresh <- as.numeric(sapply(tmp,getElement,4))
    w.beta.bar <- paste( sapply(tmp,getElement,1),sapply(tmp,getElement,2),sapply(tmp,getElement,3), sep='_' )
    u.w.thresh <- unique(w.thresh)
    p.clump <- sapply(sapply(beta.bar1,getElement,'p.value'),getElement,1)
    beta.bar11 <- as.data.frame(matrix( nrow=0, ncol=4 ))

    for( i in 1:length(beta.bar1) ){
        j <- 1
        beta.bar <- as.data.frame(matrix( nrow=nrow(beta.bar1[[i]]), ncol=4 ))
        beta.bar[,1:3] <- beta.bar1[[i]][,c('snp','effect.allele','ref.allele')]
        beta.bar[,4] <- 0
        while( j<=length(u.w.thresh) & p.clump[i] < u.w.thresh[j] ){
            ptr.weight <- which( w.thresh==u.w.thresh[j] )
            ptr.beta.bar <- match( w.beta.bar[ptr.weight], names(beta.bar1[[i]]) )
            beta.bar[,4] <- beta.bar[,4] +
                as.matrix(beta.bar1[[i]][,ptr.beta.bar]) %*% w.ridge11[ptr.weight]
            j <- j+1
        }
        if( j>1 ){
            beta.bar11 <- rbind( beta.bar11, beta.bar )
        }
    }

    if( !is.null(opt$pred2) ){
        i.snps <- intersect(beta.bar11$V1,beta.bar22$V1)
        ptr1 <- match( i.snps, beta.bar11$V1 )
        ptr2 <- match( i.snps, beta.bar22$V1 )
        beta.bar11[ ptr1, 4] <- beta.bar11[ ptr1, 4] + beta.bar22[ ptr2, 4]
        beta.bar.chr <- rbind( beta.bar11, beta.bar22[-ptr2,] )
        beta.bar.genome <- rbind( beta.bar.genome, beta.bar.chr )
    }else{
        beta.bar.genome <- rbind( beta.bar.genome, beta.bar11 )
    }
}
colnames(beta.bar.genome) <- c('snp','effect.allele','ref.allele','effect')
write.table( beta.bar.genome,wtfile,col.names=TRUE, row.names=FALSE, quote=FALSE )

if( !is.null(warnings()) ){
    print(warnings())
}

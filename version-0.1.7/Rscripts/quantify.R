library(data.table)
library(boot)
library(glmnet)
library(doMC)
library("optparse")
#source('~/BridgePRS/src/functions.R')

var.explained <- function(data,ptr){
    ptr.X <- grep('X',colnames(data))
    ptr.PRS <- grep('PRS',colnames(data))
    fit0 <- summary(lm( data$y ~ 0 + as.matrix(data[,ptr.X]), subset=ptr ))
    fit1 <- summary(lm( data$y ~ 0 + as.matrix(data[,ptr.X]) + as.matrix(data[,ptr.PRS]), subset=ptr ))
    R2 <- 1 - (1-fit1$adj.r.squared) / (1-fit0$adj.r.squared)
    return(R2)
}

option_list = list(
    make_option(c("--fpath"), type="character", default=NULL,help="", metavar="character"),
    make_option(c("--pred.dir1"), type="character", default=NULL,help="", metavar="character"),
    make_option(c("--pred.dir2"), type="character", default=NULL,help="", metavar="character"),
    make_option(c("--pred1"), type="character", default=NULL,help="", metavar="character"),
    make_option(c("--eval1"), type="character", default=NULL,help="", metavar="character"),
    make_option(c("--pred2"), type="character", default=NULL,help="", metavar="character"),
    make_option(c("--outdir"), type="character", default=NULL,help="", metavar="character"),
    make_option(c("--pop2"), type="character", default=NULL,help="", metavar="character"),
    make_option(c("--valid.data"), type="character", default=NULL,help="", metavar="character"),
    make_option(c("--test.data"), type="character", default=NULL,help="", metavar="character"),
    make_option(c("--cov.names"), type="character", default="000",help="", metavar="character"),
    make_option(c("--pheno.name"), type="character", default=NULL,help="", metavar="character"),
    make_option(c("--ids.col"), type="logical", default=TRUE,help="", metavar="logical"),
    make_option(c("--n.cores"), type="numeric", default=1,help="Number of processors for mclapply to use", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print(opt)

source(opt$fpath) 
#source('~/BridgePRS/src/functions.R')


tmp <- t(data.frame(opt))
rownames(tmp) <- names(opt)
#### Prediction ####
pTest <-     paste0(opt$pred1,'_all_preds_test.dat')
pValid <-    paste0(opt$pred1,'_all_preds_valid.dat')
outPrefix <- paste0(opt$outdir,"/",opt$pop2,"_quantify") 


logfile <- paste0(outPrefix,"_weighted_combined_var_explained.log")
write.table(tmp,file=logfile,quote=FALSE,col.names=FALSE)

if( opt$cov.names!="000" ){
    cov.names <- strsplit( opt$cov.names, ',' )[[1]]
}
registerDoMC(cores = opt$n.cores)


pred1 <- fread(pTest, data.table=FALSE) 
target <- fread(opt$test.data,data.table=FALSE)
if( opt$ids.col ){
    target <- target[match(pred1$id,target$IID),]
    pred1 <- pred1[,-1]
}



v.target <- var(target[,opt$pheno.name])
nfolds <- ifelse( nrow(target)<2000, nrow(target), 50 )
write( paste("Test data of", nrow(target), "samples"), file=logfile )
fit.ridge1 <- cv.glmnet( y=target[,opt$pheno.name], x=as.matrix(pred1),family="gaussian",alpha=0, parallel=TRUE, nfolds=nfolds, grouped=FALSE )
w.ridge1 <- getGlmnetFit( fit.ridge1, as.matrix(pred1), s='lambda.min', sparse=FALSE )[-1]



#### Validation ####
if( opt$valid.data!=0 ){
    pred1 <- fread(pValid, data.table=FALSE) 
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
    data <- data.frame( target[  ,opt$pheno.name], prs.ridge1, covs )
    my_data <- data.frame( target$IID, target[,opt$pheno.name], prs.ridge1)
    colnames(my_data) <- c('---','pheno','prs') 
    write.table( my_data, paste0(outPrefix,"_prs_predictions.out"), row.names = FALSE) 
    

    colnames(data)[1:2] <- c('y','PRS')
    b <- boot(data,var.explained,stype="i",R=10000,parallel='multicore',ncpus=opt$n.cores)
    ci <- boot.ci(b,type='norm')
    ptr.min <- which(fit.ridge1$lambda==fit.ridge1$lambda.min)
    VE.ridge1 <- c( b$t0,ci$normal[-1], fit.ridge1$cvm[ptr.min], fit.ridge1$cvsd[ptr.min] )

    if( is.null(opt$pred2) ){
        out <- rbind( VE.ridge1 )
        colnames(out) <- c('Est','2.5%','97.5%','cv.dev','cv.dev.sd')
        rownames(out) <- c('Ridge')
        print(out)
        write.csv( out, paste0(outPrefix,"_weighted_combined_var_explained.txt"), row.names=TRUE )
    }
}



tmp <- strsplit( names(w.ridge1), "_" )
lambda <- as.numeric(sapply(tmp,getElement,2))
alpha <- as.numeric(sapply(tmp,getElement,3))
lambda.weights2 <- tapply( w.ridge1, lambda, sum )
alpha.weights2 <- tapply( w.ridge1, alpha, sum )
if( is.null(opt$pred2) ){
    #write.table( alpha.weights2,paste0(opt$outdir,"/",opt$pop2,"_alpha_weights.dat"), row.names=TRUE )
    #write.table( lambda.weights2,paste0(opt$outdir,"/",opt$pop2,"_lambda_weights.dat"), row.names=TRUE )
    write.table( alpha.weights2, paste0(outPrefix,"_alpha_weights.dat"), row.names=TRUE )
    write.table( lambda.weights2,paste0(outPrefix,"_lambda_weights.dat"), row.names=TRUE )
}


w.ridge11 <- w.ridge1
beta.bar.genome <- as.data.frame(matrix( nrow=0, ncol=4 ))
for( chr in 1:22 ){
    
        
    cModel <- paste0(opt$eval1,'_beta_bar_chr',chr,'.txt.gz') 
    beta.bar1 <- fread(cModel, data.table=FALSE) 
    beta.bar1 <- beta.list(beta.bar1)

    tmp <- strsplit( names(w.ridge11), '_' )
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
    beta.bar.genome <- rbind( beta.bar.genome, beta.bar11 )


}
colnames(beta.bar.genome) <- c('snp','effect.allele','ref.allele','effect')
write.table( beta.bar.genome, paste0(outPrefix,"_weighted_combined_snp_weights.dat"), col.names=TRUE, row.names=FALSE, quote=FALSE )

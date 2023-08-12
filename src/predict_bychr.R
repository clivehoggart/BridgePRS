library(BEDMatrix)
library(MASS)
library(parallel)
library(data.table)
library("optparse")
source('~/BridgePRS/src/functions.R')
options(stringsAsFactors=FALSE)

option_list = list(
    make_option(c("-b", "--beta.stem"), type="character",
                help="Beta estimate file", metavar="character"),
    make_option(c("-o", "--out.file"), type="character",
                help="Output file stem", metavar="character"),
    make_option(c("--p.thresh"), type="character",
                default="1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8",
                help="P-value threshold to select clumps",
                metavar="character"),
    make_option(c("--non.overlapping"), type="logical", default=FALSE,
                help="Flag to set non-overlapping p-value thresholds",
                metavar="character"),
    make_option(c("--all.preds"), type="logical", default=TRUE,
                help="Flag to write all predictions", metavar="character"),
    make_option(c("--bfile"), type="character",
                help="Plink file to extract test & train genotype data",
                metavar="character"),
    make_option(c("--test.data"), type="character",
                help="File of test data and ids to use from bfile",
                metavar="character"),
    make_option(c("--valid.data"), type="character",
                help="File of validation data and ids to use from bfile",
                metavar="character"),
    make_option(c("--n.cores"), type="numeric", default=1,
                help="Number of processors for mclapply to use", metavar="character"),
    make_option(c("--cov.names"), type="character", default="000",
                help="Names of covariates to include", metavar="character"),
    make_option(c("--pheno.name"), type="character",
                help="Name of phenotype", metavar="character"),
    make_option(c("--ranking"), type="character", default="pv",
                help="Method to rank loci", metavar="character"),
    make_option(c("--by.chr"), type="numeric", default=1,
                help="Logical indicating if bed files are split by chr",
                metavar="character"),
    make_option(c("--strand.check"), type="numeric", default=0,
                help="Keep only non-ambiguous SNPs", metavar="numeric"),
    make_option(c("--binary"), type="numeric", default=0,
                help="Indicator for binary outcome", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

print(opt)

if( opt$binary==1 ){
    family <- "binomial"
}else{
    family <- "gaussian"
}

logfile <- paste0(opt$out.file,".log")
tmp <- t(data.frame(opt))
rownames(tmp) <- names(opt)
write.table(tmp,file=logfile,quote=FALSE,col.names=FALSE)

if( opt$cov.names!="000" ){
    cov.names <- strsplit( opt$cov.names, ',' )[[1]]
}
test.data <- fread(opt$test.data,data.table=FALSE)
ids.use <- as.character( test.data$IID )

if( opt$valid.data!=0 ){
    valid.data <- fread(opt$valid.data,data.table=FALSE)
    ids.use <- as.character(c( ids.use, valid.data$IID ))
}

tmp <- as.numeric(strsplit( opt$p.thresh, ',' )[[1]])
p.thresh <- matrix(ncol=2,nrow=length(tmp))
p.thresh[,1] <- tmp
if( opt$non.overlapping ){
    p.thresh[,2] <- c(tmp[-1],-1)
}else{
    p.thresh[,2] <- -1
}

if( opt$by.chr==0 ){
    ptr.bed <- BEDMatrix( opt$bfile, simple_names=TRUE )
    bim <- fread( paste(opt$bfile,'.bim',sep='' ) )
    ids.use <- intersect( ids.use, attributes(ptr.bed)[[3]][[1]] )
}
for( chr in 1:22 ){
    if( opt$by.chr==1 ){
        ptr.bed <- BEDMatrix( paste(opt$bfile,chr,sep=''), simple_names=TRUE )
        bim <- fread( paste(opt$bfile,chr,'.bim',sep='' ) )
        ids.use <- intersect( ids.use, attributes(ptr.bed)[[3]][[1]] )
    }
#    beta.bar <- dget(paste(opt$beta.file,chr,'.Robj',sep=''))
    beta.bar <- fread(paste0(opt$beta.stem,"_beta_bar_chr",chr,".txt.gz"),data.table=FALSE)
    clump.ids <- unique(beta.bar$clump.id)
    beta.bar <- beta.list(beta.bar)
    names(beta.bar) <- clump.ids
#Use same thresholds for all chromosomes miroring p-value distribution
    if( opt$ranking=='f.stat' | opt$ranking=='thinned.f.stat' ){
        infile <- paste0(opt$beta.stem,"_KLdist_chr",chr,".txt.gz")
        kl.metric <- fread(infile, data.table=FALSE)[,-1]
        if( chr==1 ){
            p.clump <- sapply(sapply(beta.bar,getElement,'p.value'),getElement,1)
            F <- ecdf(p.clump)
            q <- F(10^(-(1:8)))
            kl.thresh <- apply( kl.metric, 2, quantile, 1-q, na.rm=TRUE )
            n.thresh <- nrow(kl.thresh)
            nme.thresh <- rownames(kl.thresh)
        }
    }
    if( opt$ranking=='pv' | opt$ranking=='pv.minP' |
        opt$ranking=='pv.ftest' | opt$ranking=='thinned.pv.ftest' ){
        kl.metric <- NULL
        kl.thresh <- NULL
        n.thresh <- nrow(p.thresh)
        nme.thresh <- p.thresh[,1]
    }

    pred.chr <- get.pred.genome( beta.bar,
                                p.thresh=p.thresh, X.bed=ptr.bed,
                                bim=bim, ids.use=ids.use, b.size=100,
                                n.cores=as.numeric(opt$n.cores),
                                kl.metric, kl.thresh, opt$ranking, opt$strand.check )
    if( chr==1 ){
        pred.genome <- pred.chr
    }
    if( chr>1 ){
        for( j in 1:n.thresh ){
            pred.genome[[j]] <- pred.genome[[j]] + pred.chr[[j]]
        }
    }
}

if( opt$valid.data!=0 ){
    VE.valid <- matrix(ncol=n.thresh,nrow=ncol(pred.genome[[1]]))
    rownames(VE.valid) <- colnames(pred.genome[[1]])
    colnames(VE.valid) <- nme.thresh

    valid.ids <- intersect( valid.data$IID, ids.use )
    valid.data <- valid.data[match(valid.ids, valid.data$IID ),]
    if( opt$cov.names!="000" ){
        valid.covs <- cbind( 1, as.matrix(valid.data[,cov.names]) )
    }else{
        valid.covs <- rep(1,nrow(valid.data))
    }
    ptr.valid <- match( valid.data$IID, rownames(pred.genome[[1]]) )
    all.preds.valid <- as.data.table(valid.ids)
    colnames(all.preds.valid)[1] <- 'ids'

    if( family=="gaussian" )
        fit.valid0 <- summary(lm( valid.data[,opt$pheno.name] ~ 0 + valid.covs ))
    else
        fit.valid0 <- summary(glm( valid.data[,opt$pheno.name] ~ 0 + valid.covs, family=family ))
    }

VE.test <- matrix(ncol=n.thresh,nrow=ncol(pred.genome[[1]]))
rownames(VE.test) <- colnames(pred.genome[[1]])
colnames(VE.test) <- nme.thresh

test.ids <- intersect( test.data$IID, ids.use )
test.data <- test.data[match(test.ids, test.data$IID ),]
if( opt$cov.names!="000" ){
    test.covs <- cbind( 1, as.matrix(test.data[,cov.names]) )
}else{
    test.covs <- rep(1,nrow(test.data))
}
ptr.test <- match( test.data$IID, rownames(pred.genome[[1]]) )
all.preds.test <- as.data.table(test.ids)
colnames(all.preds.test)[1] <- 'ids'

if( family=="gaussian" )
    fit.test0 <- summary(lm( test.data[,opt$pheno.name] ~ 0 + test.covs ))
else
    fit.test0 <- summary(glm( test.data[,opt$pheno.name] ~ 0 + test.covs, family=family ))

for( k in 1:n.thresh ){
    for( j in 1:ncol(pred.genome[[k]]) ){
        if( family=="gaussian" ){
            fit.test1 <- summary(lm( test.data[,opt$pheno.name] ~
                                         0 + pred.genome[[k]][ptr.test,j] + test.covs ))
            VE.test[j,k] <- 1 - (1-fit.test1$adj.r.squared) / (1-fit.test0$adj.r.squared)
            if( opt$valid.data!=0 ){
                fit.valid1 <- summary(lm( valid.data[,opt$pheno.name] ~
                                              0 + pred.genome[[k]][ptr.valid,j] + valid.covs ))
                VE.valid[j,k] <- 1 - (1-fit.valid1$adj.r.squared) / (1-fit.valid0$adj.r.squared)
            }
        }else{
            fit.test1 <- summary(glm( test.data[,opt$pheno.name] ~
                                          0 + pred.genome[[k]][ptr.test,j] + test.covs,
                                     family=family ))
            n <- length(fit.test0$y)
            d_null <- -2 * logLik(fit.test0)[1]
            d_full <- -2 * logLik(fit.test1)[1]
            VE.test[j,k] <- (1 - exp((d_full - d_null) / n)) / (1 - exp(-d_null / n))
            if( opt$valid.data!=0 ){
                fit.valid1 <- summary(hlm( valid.data[,opt$pheno.name] ~
                                               0 + pred.genome[[k]][ptr.valid,j] + valid.covs,
                                          family=family ))
                n <- length(fit.valid0$y)
                d_null <- -2 * logLik(fit.valid0)[1]
                d_full <- -2 * logLik(fit.valid1)[1]
                VE.valid[j,k] <- (1 - exp((d_full - d_null) / n)) / (1 - exp(-d_null / n))
            }
        }
    }
    if( opt$all.preds ){
        if( length(colnames(pred.genome[[k]]))>0 ){
            tmp <- pred.genome[[k]]
            colnames(tmp) <- paste( colnames(tmp), nme.thresh[k], sep='_' )
            all.preds.test <- cbind( all.preds.test, tmp[ptr.test,] )
            if( opt$valid.data!=0 ){
                all.preds.valid <- cbind( all.preds.valid, tmp[ptr.valid,] )
            }
        }
    }
}
ptr.col <- order(apply(VE.test,2,max),decreasing=TRUE)[1]
s <- order(VE.test[,ptr.col],decreasing=TRUE)

if( opt$all.preds ){
    fwrite( all.preds.test, paste( opt$out.file, '_all_preds_test.dat',sep='' ))
    if( opt$valid.data!=0 ){
        fwrite( all.preds.valid, paste( opt$out.file, '_all_preds_valid.dat',sep=''))
    }
}

if( opt$valid.data!=0 ){
    write.table( VE.valid[s,], paste( opt$out.file, '_var_exp_valid.dat', sep='' ), quote=FALSE )
    out <- data.frame( rownames(pred.genome[[1]])[ptr.valid], pred.genome[[ptr.col]][ptr.valid,s[1]] )
    colnames(out) <- c('IID','PRS')
    fwrite( out, paste( opt$out.file, '_best_pred_valid.dat', sep='' ), quote=FALSE )
}

write.table( VE.test[s,], paste( opt$out.file, '_var_exp_test.dat', sep='' ), quote=FALSE )
out <- data.frame( rownames(pred.genome[[1]])[ptr.test], pred.genome[[ptr.col]][ptr.test,s[1]] )
colnames(out) <- c('IID','PRS')
fwrite( out, paste( opt$out.file, '_best_pred_test.dat', sep='' ), quote=FALSE )

tmp <-  strsplit( colnames(pred.genome[[1]])[s[1]], '_' )[[1]]
if( length(tmp)==3 ){
    lambda.opt <- as.numeric(tmp[[2]])
    S.opt <- as.numeric(tmp[[3]])
    p.opt <- nme.thresh[ptr.col]
    write.table( data.frame(S.opt,lambda.opt,p.opt),
                paste( opt$out.file, '_best_model_params.dat', sep='' ),
                row.names=FALSE, quote=FALSE )
}
if( length(tmp)==2 ){
    tau.opt <- as.numeric(tmp[[2]])
    p.opt <- nme.thresh[ptr.col]
    write.table( data.frame(tau.opt,p.opt),
                paste( opt$out.file, '_best_model_params.dat', sep='' ),
                row.names=FALSE, quote=FALSE )
}
#prsice <- fread('/sc/arion/work/hoggac01/platelets/prsice.best')
#ptr <- match(test.data$IID,prsice$IID)
#fit.test1 <- summary(lm( test.data[,opt$pheno.name] ~ prsice$PRS[ptr] +
#                              as.matrix(test.data[,cov.names]) ))

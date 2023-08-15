library(BEDMatrix)
library(MASS)
library(parallel)
library(data.table)
library("optparse")
source('~/BridgePRS/bin/functions.R')
options(stringsAsFactors=FALSE)

opt <- list()
opt$beta.stem <- "/sc/arion/projects/psychgen/projects/prs/cross_population_prs_development/quick_ridge/results/BBJ/test/50/models/EAS_stage1"
opt$out.file <- "/sc/arion/projects/psychgen/projects/prs/cross_population_prs_development/quick_ridge/results/BBJ/test/50/EAS_stage1"
opt$p.thresh <- "1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8"
opt$non.overlapping <- TRUE
opt$all.preds <- TRUE
opt$bfile <- "/sc/arion/projects/psychgen/ukb/usr/clive/ukb/imputed/chr"
opt$test.data <- "/sc/arion/projects/psychgen/projects/prs/cross_population_prs_development/quick_ridge/hapmap/phenotype/50/EAS-50-test.dat"
opt$valid.data <- "/sc/arion/projects/psychgen/projects/prs/cross_population_prs_development/quick_ridge/hapmap/phenotype/50/EAS-50-valid.dat"
opt$n.cores <- 30
opt$cov.names <- "Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21,PC22,PC23,PC24,PC25,PC26,PC27,PC28,PC29,PC30,PC31,PC32,PC33,PC34,PC35,PC36,PC37,PC38,PC39,PC40"
opt$pheno.name <- "Phenotype"
opt$by.chr <- 1
opt$strand.check <- 1
opt$help <- FALSE

if( opt$cov.names!="000" ){
    cov.names <- strsplit( opt$cov.names, ',' )[[1]]
}
test.data <- fread(opt$test.data,data.table=FALSE)
ids.use <- as.character( test.data$IID )
if( opt$cov.names!="000" ){
    test.covs <- cbind( 1, as.matrix(test.data[,cov.names]) )
}else{
    test.covs <- rep(1,nrow(test.data))
}
if( opt$valid.data!=0 ){
    valid.data <- fread(opt$valid.data,data.table=FALSE)
    ids.use <- as.character(c( ids.use, valid.data$IID ))
    if( opt$cov.names!="000" ){
        valid.covs <- cbind( 1, as.matrix(valid.data[,cov.names]) )
    }else{
        valid.covs <- rep(1,nrow(valid.data))
    }
}

if( opt$by.chr==0 ){
    ptr.bed <- BEDMatrix( opt$bfile, simple_names=TRUE )
    bim <- fread( paste(opt$bfile,'.bim',sep='' ) )
    ids.use <- intersect( ids.use, attributes(ptr.bed)[[3]][[1]] )
}


VE.valid <- matrix(ncol=30*22,nrow=55)
valid.ids <- intersect( valid.data$IID, ids.use )
valid.data <- valid.data[match(valid.ids, valid.data$IID ),]
fit.valid0 <- summary(lm( valid.data[,opt$pheno.name] ~ 0 + valid.covs ))

VE.test <- matrix(ncol=30*22,nrow=55)
test.ids <- intersect( test.data$IID, ids.use )
test.data <- test.data[match(test.ids, test.data$IID ),]
fit.test0 <- summary(lm( test.data[,opt$pheno.name] ~ 0 + test.covs ))

k.clump <- vector()
p.clump <- vector()
for( chr in 1:22 ){
    beta.bar <- fread(paste0(opt$beta.stem,"_beta_bar_chr",chr,".txt.gz"),data.table=FALSE)
    clump.ids <- unique(beta.bar$clump.id)
    beta.bar <- beta.list(beta.bar)
    names(beta.bar) <- clump.ids
    ptr.beta.use <- grep( 'beta.bar', colnames(beta.bar[[1]]) )
    if( opt$by.chr==1 ){
        ptr.bed <- BEDMatrix( paste(opt$bfile,chr,sep=''), simple_names=TRUE )
        bim <- fread( paste(opt$bfile,chr,'.bim',sep='' ) )
        ids.use <- intersect( ids.use, attributes(ptr.bed)[[3]][[1]] )
    }
    for( i in 1:30 ){
        ii <-  i + (chr-1)*30
        k.clump[ii] <- nrow(beta.bar[[i]])
        p.clump[ii] <- beta.bar[[i]]$p.value[1]
        pred.clump <- get.pred.clump( beta.bar[[i]], clump.id=names(beta.bar)[i],
                                     ptr.beta.use=ptr.beta.use, by.chr=0,
                                     X.bed=ptr.bed, bim=bim, ids.use=ids.use,
                                     opt$strand.check )

        rownames(VE.valid) <- colnames(pred.clump)
        ptr.valid <- match( valid.data$IID, rownames(pred.clump) )
        rownames(VE.test) <- colnames(pred.clump)
        ptr.test <- match( test.data$IID, rownames(pred.clump) )

        for( j in 1:ncol(pred.clump) ){
            fit.test1 <- summary(lm( test.data[,opt$pheno.name] ~
                                         0 + pred.clump[ptr.test,j] + test.covs ))
            VE.test[j,ii] <- 1 - (1-fit.test1$adj.r.squared) / (1-fit.test0$adj.r.squared)
            if( opt$valid.data!=0 ){
                fit.valid1 <- summary(lm( valid.data[,opt$pheno.name] ~
                                              0 + pred.clump[ptr.valid,j] + valid.covs ))
                VE.valid[j,ii] <- 1 - (1-fit.valid1$adj.r.squared) / (1-fit.valid0$adj.r.squared)
            }
        }
    }
}
VE.both <- VE.test + VE.valid
tmp <- strsplit( rownames(VE.test), '_' )
lambda <- sapply( tmp, getElement, 2 )
s0 <- apply(VE.both,2,order,decreasing=TRUE)
s1 <- apply(VE.test,2,order,decreasing=TRUE)
s2 <- apply(VE.valid,2,order,decreasing=TRUE)
VE.mx.test <- apply( VE.test, 2, max )
VE.mx.valid <- apply( VE.valid, 2, max )
ptr <- which( VE.mx.test>0 & VE.mx.valid>0 )
for( i in ptr ){
    print(c(i,k.clump[i]))
    print(VE.both[s0[1:4,i],i])
    print(VE.test[s1[1:4,i],i])
    print(VE.valid[s2[1:4,i],i])
    print(VE.test[s2[1:4,i],i])
    print(VE.valid[s1[1:4,i],i])
}
q.p.clump <- assign.to.quants( p.clump, 5 )
q.k.clump <- assign.to.quants( k.clump, 5 )

sum.p <- matrix(nrow=5,ncol=55)
sum.k <- matrix(nrow=5,ncol=55)
for( i in 1:5 ){
    sum.p[i,] <- apply(VE.both[,q.k.clump==(i-1)],1,sum)
    sum.k[i,] <- apply(VE.both[,q.p.clump==(i-1)],1,sum)
}
colnames(sum.p) <- rownames(VE.both)
colnames(sum.k) <- rownames(VE.both)
s.p <- apply(sum.p,1,order,decreasing=TRUE)
s.k <- apply(sum.k,1,order,decreasing=TRUE)
for( i in 1:5 ){
    print( sum.p[i,s.p[1:4,i]] )
}
for( i in 1:5 ){
    print( sum.k[i,s.k[1:4,i]] )
}
tapply(p.clump,q.p.clump,range)
tapply(k.clump,q.k.clump,range)
q.p.clump <- assign.to.quants( p.clump, 5 )
q.k.clump <- assign.to.quants( k.clump, 5 )

sum.p <- matrix(nrow=5,ncol=55)
sum.k <- matrix(nrow=5,ncol=55)
for( i in 1:5 ){
    sum.p[i,] <- apply(VE.both[,q.k.clump==(i-1)],1,sum)
    sum.k[i,] <- apply(VE.both[,q.p.clump==(i-1)],1,sum)
}
colnames(sum.p) <- rownames(VE.both)
colnames(sum.k) <- rownames(VE.both)
s.p <- apply(sum.p,1,order,decreasing=TRUE)
s.k <- apply(sum.k,1,order,decreasing=TRUE)
for( i in 1:5 ){
#    print( sum.p[i,s.p[1:4,i]] )
    print( 100*max(apply(VE.both[,q.p.clump==(i-1)],1,mean)) )
}
for( i in 1:5 ){
#    print( sum.k[i,s.k[1:4,i]] )
    print( 100*max(apply(VE.both[,q.k.clump==(i-1)],1,mean)) )
}
tapply(p.clump,q.p.clump,median)
tapply(k.clump,q.k.clump,median)
table(q.p.clump,q.k.clump)

sum.k <- matrix(nrow=5,ncol=55)
for( i in 1:5 ){
    sum.p[i,] <- apply(VE.both[,q.k.clump==(i-1)],1,sum)
    sum.k[i,] <- apply(VE.both[,q.p.clump==(i-1)],1,sum)
}

for( i in 0:4 ){
    tmp0 <- vector()
    tmp <- vector()
    for( j in 0:4 ){
        ptr <- which( q.p.clump==i & q.k.clump==j )
        vec <- apply(VE.both[,ptr],1,mean)
        tmp0[j+1] <- 100*max(vec)
        s <- order( vec, decreasing=TRUE )
        tmp[j+1] <- 100*vec[s[1]]
        tmp0[j+1] <- lambda[s[1]]

    }
    print(tmp0)
    print(tmp)
}

ptr <- which( q.p.clump==0 & q.k.clump==4 )
vec <- apply(VE.both[,ptr],1,mean)
s <- order( vec, decreasing=TRUE )
vec[s]
ptr <- which( q.p.clump==0 & q.k.clump==3 )
vec <- apply(VE.both[,ptr],1,mean)
s <- order( vec, decreasing=TRUE )
vec[s]

regularize_ld_matrix <- function(ld.mat, target_min_ev = NULL, verbose = TRUE) {
  stopifnot(is.matrix(ld.mat), isSymmetric(ld.mat))

  # Compute eigenvalues once
  ev <- eigen(ld.mat, symmetric = TRUE, only.values = TRUE)$values
  min_ev <- min(ev)

  # Default threshold based on matrix scale
  if (is.null(target_min_ev)) {
    target_min_ev <- max(1e-4, 0.01 * mean(ev))
  }

  if (min_ev > target_min_ev) {
    if (verbose) cat("Matrix already positive definite. No regularization needed.\n")
    return(ld.mat)
  }

  # Compute required lambda
  lambda <- target_min_ev - min_ev
  ld.mat.reg <- ld.mat + lambda * diag(diag(ld.mat))

  # Check positive definiteness via Cholesky
  success <- tryCatch({
    chol(ld.mat.reg)
    TRUE
  }, error = function(e) FALSE)

  if (!success) {
    if (verbose) cat("Increasing lambda to ensure Cholesky success...\n")
    factor <- 10
    count <- 0
    while (!success && count < 5) {
      lambda <- lambda * factor
      ld.mat.reg <- ld.mat + lambda * diag(diag(ld.mat))
      success <- tryCatch({
        chol(ld.mat.reg)
        TRUE
      }, error = function(e) FALSE)
      count <- count + 1
    }

    if (!success) stop("Failed to regularize LD matrix to positive definiteness.")
  }

  if (verbose) cat("Final lambda used:", format(lambda, scientific = TRUE), "\n")

  return(ld.mat.reg)
}

standardise <- function(X){
    m <- apply(X,2,mean,na.rm=TRUE)
    X <- t(t(X)-m)
    s <- apply(X,2,sd,na.rm=TRUE)
    X <- t(t(X)/s)
    return(X)
}

allele.check <- function( effect.allele1, ref.allele1,
                         effect.allele2, ref.allele2, strand.check ){
    swtch <- ifelse( effect.allele1==effect.allele2 &
                     ref.allele1==ref.allele2, 1, 0 )
    swtch <- ifelse( effect.allele1==ref.allele2 &
                     ref.allele1==effect.allele2, -1, swtch )
    ptr.miss <- which( swtch==0 )
    if( strand.check & length(ptr.miss)>0 ){
        effect.allele22 <- alt.strand( effect.allele2 )
        ref.allele22 <- alt.strand( ref.allele2 )

        swtch[ptr.miss] <- ifelse( effect.allele1[ptr.miss]==effect.allele22[ptr.miss] &
                               ref.allele1[ptr.miss]==ref.allele22[ptr.miss],
                               1, swtch[ptr.miss] )

        swtch[ptr.miss] <- ifelse( effect.allele1[ptr.miss]==ref.allele22[ptr.miss] &
                                   ref.allele1[ptr.miss]==effect.allele22[ptr.miss],
                                  -1, swtch[ptr.miss] )
    }
    return(swtch)
}

genos.cor <- function(genos,maf.thresh=0.005){
    m <- apply(genos,2,mean,na.rm=TRUE)
    maf <- ifelse( m<1, m, 2-m ) / 2
    genos <- t(t(genos)-m)
    genos <- genos[,which(maf>maf.thresh),drop=FALSE]
    s <- apply(genos,2,sd,na.rm=TRUE)
    genos <- t(t(genos)/s)
    genos <- ifelse( is.na(genos), 0, genos )
    M <- cor(genos,use="pairwise.complete.obs")
    return(M)
}

eff.gao <- function(M=NULL,genos=NULL,X=NULL){
    if( !is.null(genos) ){
        M <- genos.cor(genos)
    }
    if( !is.null(X) ){
        M <- cor(X,use="pairwise.complete.obs")
    }
    T <- eigen(M,only.values=TRUE)
    T <- T$values
    TT <- cumsum(T)/sum(T)
    ef <- min(which(TT>0.995))
    return(ef)
}

beta.tail.approx <- function( x1, alpha, beta ){
# Trapezium approximation
    inc <- 1/(4*beta)
    x <- seq(from=x1,by=inc,length=100)
    logf <- (alpha-1)*log(x) + (beta-1)*log(1-x) - lgamma(alpha) - lgamma(beta) + lgamma(alpha+beta)
    mx.logf <- max(logf)
    f <- exp(logf-mx.logf)
    height <- (f[-100] + f[-1]) / 2
    area <- log(sum(height*inc)) + mx.logf
    return(area)
}

alt.strand <- function(allele){
    allele1 <- ifelse( allele=='A', 'T', allele )
    allele1 <- ifelse( allele=='C', 'G', allele1 )
    allele1 <- ifelse( allele=='G', 'C', allele1 )
    allele1 <- ifelse( allele=='T', 'A', allele1 )
    return(allele1)
}

getGlmnetFit <- function( fit, X, s="lambda.min", sparse=TRUE ){
    tmp2 <- coef( fit, X, s=s )
    beta <- as.vector(tmp2)
    names(beta) <- rownames(tmp2)
    if( sparse ){
        ptr <- which( tmp2!=0 )
        beta <- beta[ptr]
    }
    return(beta)
}

ld.shrink <- function( snps, bim, recomb, m, Ne ){
    posn <- bim$V4[match( snps, bim$V2 )]

    ptr <- match( posn, recomb$Position )

    ptr.begin <- which( recomb$Position < min(posn,na.rm=TRUE) )
    ptr.begin <- ifelse( length(ptr.begin)>1, max(ptr.begin), 1 )
    ptr.end <- which( recomb$Position > max(posn,na.rm=TRUE) )
    ptr.end <- ifelse( length(ptr.end)>1, min(ptr.end), nrow(recomb) )

    map <- c( recomb$Map[ ptr.begin ], recomb$Map[ptr], recomb$Map[ ptr.end ] )
    start.posn <- min(c( recomb$Position[ptr.begin], posn ))
    end.posn <- max(c( recomb$Position[ptr.end], posn ))

    posn <- c( start.posn, posn, end.posn )

    s <- order(posn)
    ss <- order(s)
    cbind(posn,map)[s,]
    map <- map[s]
    posn <- posn[s]
    ptr0 <- which(is.na(map))
    ptr1 <- which(!is.na(map))

    for( i in 1:(length(ptr1)-1) ){
        ptr.fill <- which( ptr1[i]< ptr0&ptr0 < ptr1[i+1] )
        if( length(ptr.fill)>0 ){
            a <- map[ptr1[i]]
            b <- (map[ptr1[i+1]] - map[ptr1[i]]) / (posn[ptr1[i+1]] - posn[ptr1[i]])
            map[ptr0[ptr.fill]] <- a + b*(posn[ptr0[ptr.fill]] - posn[ptr1[i]] )
        }
    }
    map <- map[ss]
    posn <- posn[ss]
    cbind(posn,map)[s,]

    map <- map[-c(1,length(map))]
    ld.shrink.factor <- diag(x=1,nrow=length(map),ncol=length(map))
    for( i in 1:(length(map)-1) ){
        d <- abs(map[-(1:i)] - map[i])
        dd <- exp(-2*Ne*d/m)
        ld.shrink.factor[i,-(1:i)] <- dd
        ld.shrink.factor[-(1:i),i] <- dd
    }
    return( ld.shrink.factor )
}

beta.list <- function( beta ){
    clump.id <- unique(beta$clump.id)
    newbeta <- list(beta)
    for( i in 1:length(clump.id) ){
        ptr <- which( beta$clump.id==clump.id[i] )
        newbeta[[i]] <- beta[ptr,]
    }
    return(newbeta)
}

kl.dist <- function( mu0, mu1, lambda0, lambda1, n ){
    k <- length(mu0)
    KL.01 <- ( log(det(lambda1)/det(lambda0))
        + sum(diag((solve(lambda1)%*%lambda0)))
        + t(mu1-mu0) %*% (n*lambda0) %*% (mu1-mu0) - k ) / 2
    return(KL.01)
}

noncentral.ridge.fit <- function( beta.data, LD, af,
                                 beta.tilde1, lambda0, lambda1, w.prior,
                                 n, precision=FALSE, ranking, sigma2 ){
    s2 <- 2*af*(1-af)
#    s <- as.vector(sqrt(s2))
    k <- length(beta.data)
    ret <- vector("list", 3)

# "Ref" prior
#    beta0 <- rep(0,k)
#    lambda0 <- diag(1,nrow=k)
# Data driven posterior
#    lambda02 <- diag(s,nrow=k) %*% LD %*% diag(s,nrow=k) + lambda0
#    beta.tilde02 <- solve(lambda0) %*% as.matrix( beta.data * s*s )

#    XtX <- diag(s,nrow=k) %*% LD %*% diag(s,nrow=k)
    lambda2 <- LD + lambda1*w.prior

    beta.tilde2 <- solve(lambda2) %*% ( w.prior*lambda1%*%beta.tilde1 +
                                        as.matrix(beta.data * s2) )

# Is data driven posterior closer to informative prior than "ref" prior ?
#    d.post.prior <- kl.dist( beta.tilde02, beta.tilde1, lambda02, w.prior*lambda1, n )
#    d.post.ref <- kl.dist( beta.tilde02, beta0, lambda02, lambda0, n )
#    ret[[3]] <- d.post.ref / d.post.prior
#    ret[[3]] <- d.post.ref - d.post.prior

#    kl.02 <- kl.dist( beta.tilde2, beta0, lambda2, lambda0, n )
#    kl.12 <- kl.dist( beta.tilde2, beta.tilde1, lambda2, lambda1, n )
#    ret[[3]] <- kl.02 - kl.12

    if( ranking=="f.stat" ){
        ret[[3]] <- -Pseudo.f.test( beta.tilde2, lambda2,
                                   n.eff=n*(1+w.prior), sigma2=sigma2 )
    }

    if( ranking=="thinned.f.stat" ){
        ret[[3]] <- -Pseudo.f.test.diag( beta=beta.tilde2, lambda=lambda2,
                                        LD=ref.stats$LD,
                                        n.eff=n*(1+w.prior), sigma2=sigma2 )
    }

    ret[[1]] <- beta.tilde2
    if( precision ){
        ret[[2]] <- lambda2
    }
    return( ret )
}

ridge.fit <- function( beta.data, LD, af, l, S=1, precision=FALSE ){
    s2 <- 2*af*(1-af)
    l <- l*s2^S
#    s <- as.vector(sqrt(s2))
    k <- length(beta.data)

    beta0 <- rep(0,k)
    lambda0 <- diag(l,nrow=k)

#    lambda1 <- diag(s,nrow=k) %*% LD %*% diag(s,nrow=k) + lambda0
    lambda1 <- LD + lambda0
    beta.tilde <- solve(lambda1) %*% as.matrix( beta.data * s2 )

#    kl1 <- kl.dist( beta0, beta.tilde, lambda0, lambda1 )
#    kl2 <- kl.dist( beta.tilde, beta0, lambda1, lambda0, n )

    ret <- vector("list", 2)
    ret[[1]] <- beta.tilde
    if( precision ){
        ret[[2]] <- lambda1
    }
    return( ret )
}

f.test <- function( beta, LD, af, n, sigma2 ){
    l <- 0.01
    k <- length(beta)
    s2 <- 2*af*(1-af)

    lambda1 <- LD + diag(l,nrow=k)
    beta.hat <- solve(lambda1) %*% as.matrix( beta * s2 )
    xx <- t(beta.hat) %*% lambda1 %*% beta.hat
    stat <- (n-k) * xx / (k*sigma2)
    f.tail <- pf( stat, k, n-k, lower.tail=FALSE, log.p=FALSE )

    return( f.tail )
}

Pseudo.f.test <- function( beta, lambda, n.eff, sigma2 ){
    k <- length(beta)

    stat <- (n.eff-k) * ( t(beta) %*% lambda %*% beta ) / (k * sigma2)
#    f.tail <- pf( stat, k, n.eff-k, lower.tail=FALSE, log.p=TRUE )

    beta.stat <- as.numeric( stat*k / ( stat*k + n.eff - k ) )
    beta.tail <- pbeta( beta.stat, k/2 , (n.eff-k)/2, lower.tail=FALSE, log.p=TRUE )
    if( is.infinite(beta.tail) ){
        beta.tail <- beta.tail.approx( beta.stat, k/2 , (n.eff-k)/2 )
    }

    return( beta.tail )
}

f.test.diag <- function( beta, LD, af, n, sigma2 ){
    l <- 0.01
    k <- length(beta)
    s2 <- 2*af*(1-af)

    lambda1 <- LD + diag(l,nrow=k)
    beta.hat <- solve(lambda1) %*% as.matrix( beta * s2 )

    e <- eigen(LD,symmetric=TRUE)
    TT <- cumsum(e$values)/sum(e$values)
    k.eff <- min(which(TT>0.95))

    b <- (t(e$vector) %*% beta.hat)[1:k.eff]
    lambda <- e$values[1:k.eff]

    stat <- (n-k.eff) * sum(b * lambda * b) / (k.eff*sigma2)
    f.tail <- pf( stat, k.eff, n-k.eff, lower.tail=FALSE, log.p=FALSE )

    return( f.tail )
}

Pseudo.f.test.diag <- function( beta, lambda, LD, n.eff, sigma2 ){
    e <- eigen(LD,symmetric=TRUE)
    TT <- cumsum(e$values)/sum(e$values)
    k.eff <- min(which(TT>0.95))

    e <- eigen(lambda,symmetric=TRUE)
    b <- (t(e$vector) %*% beta)[1:k.eff]

    stat <- (n.eff-k.eff) * sum(b * e[1:k.eff] * b ) / (k.eff * sigma2)
    f.tail <- pf( stat, k.eff, n.eff-k.eff, lower.tail=FALSE, log.p=TRUE )

    return( f.tail )
}

est.ref.stats <- function( snps, ids, X.bed, bim,
                          effect.allele, ref.allele, strand.check,
                          n.eff=FALSE ){
    ids <- intersect( ld.ids, dimnames(X.bed)[[1]] )
    if( length(ids)<10 ){
        print("Fewer than 10 individuals to estimate LD")
        stop()
    }
    X <- X.bed[ ids, snps, drop=FALSE ]
    ptr.bed <- match( snps, colnames(X.bed) )
    if( length(ptr.bed)==0 ){
        print("SNPs used in clumping and ref data do not overlap")
        stop()
    }
    X <- X.bed[ ids, snps, drop=FALSE ]
    ptr.bed <- match( snps, colnames(X.bed) )
    X <- X[,match( snps, colnames(X)), drop=FALSE ]
    coded.allele <- bim$V5[ptr.bed]
    other.allele <- bim$V6[ptr.bed]

    swtch <- allele.check( coded.allele, other.allele,
                          effect.allele, ref.allele, strand.check )
#    swtch <- ifelse( coded.allele==effect.allele & other.allele==ref.allele, 1, 0 )
#    swtch <- ifelse( coded.allele==ref.allele & other.allele==effect.allele, -1, swtch )
#    ptr.miss <- which( swtch==0 )
#    if( strand.check & length(ptr.miss)>0 ){
#        coded.allele1 <- alt.strand( coded.allele )
#        other.allele1 <- alt.strand( other.allele )

#        swtch[ptr.miss] <- ifelse( coded.allele1[ptr.miss]==effect.allele[ptr.miss] &
#                                   other.allele1[ptr.miss]==ref.allele[ptr.miss],
#                                  1, 0 )

#        swtch[ptr.miss] <- ifelse( coded.allele1[ptr.miss]==ref.allele[ptr.miss] &
#                                   other.allele1[ptr.miss]==effect.allele[ptr.miss],
#                                  -1, swtch[ptr.miss] )
#    }

    # Map ref stat genotypes to effect.allele & ref.allele
    ptr <- which(swtch == -1)
    if( length(ptr)>0 ){
        X[,ptr] <- 2 - X[,ptr]
    }

    ptr <- which(swtch == 0)
    if( length(ptr)>0 ){
        X[,ptr] <- 0
    }

    af <- apply(X,2,mean,na.rm=TRUE)/2
    if( min(af)==0 ){
        print("Inclusion of monomorphic variants in clump in data used to estimate LD.")
        print("Has clumping been run independently of the Bridge package?")
    }
    ld <- cov(X,use='pairwise.complete')

    if( n.eff ){
        if( ncol(X)>1 & sum( 0.005 < af&af < 0.995 )>1 ){
            n.eff <- eff.gao(genos=X)
        }else{
            n.eff <- 1
        }
    }

    ret <- list()
    ret$af <- af
    ret$ld <- ld
    ret$n.eff <- n.eff
    return(ret)
}

sumstat.subset <- function( block.i=NULL, snp=NULL, sumstats, ld.ids,
                           X.bed, bim, n.all, n.prop, strand.check, n.folds=1, tol=1e-6 ){
    ld.shrink <- c(1e-3, 5e-3, 1e-2, 5e-2, 1e-1)
    if( !is.null(block.i) ){
        block.snps <-  unlist(strsplit( block.i$SNPS, '\\|' ))
        snps <- intersect( sumstats$SNP, block.snps )
    }
    if( !is.null(snp) ){
        snps <- intersect( sumstats$SNP, snp )
    }

    if( length(snps)>0 ){
        sumstats <- sumstats[match( snps, sumstats$SNP ),]
        ref.stats <- est.ref.stats( snps, ld.ids, X.bed, bim,
                                   sumstats$ALLELE1, sumstats$ALLELE0,
                                   strand.check, n.eff=FALSE )
        ptr.use <- which( !is.na(ref.stats$af) & ref.stats$af!=0 & ref.stats$af!=1 )
        k <- length(ptr.use)
        if( k>0 ){
            sumstats <- sumstats[ptr.use,,drop=FALSE]
            ref.stats$af <- ref.stats$af[ptr.use,drop=FALSE]
            ref.stats$ld <- ref.stats$ld[ptr.use,ptr.use,drop=FALSE]

            ld.mat <- regularize_ld_matrix( ref.stats$ld )

#            ev <- eigen(ld.mat, symmetric = TRUE, only.values = TRUE)$values
#            lambda1 <- ifelse( min(ev) < 1e-6, abs(min(ev)) + 1e-6, 0 )
#            lambda2 <- ifelse( min(ev) < 1e-6,  1e-6 - min(ev), 0 )
#            ld.mat <- ld.mat + lambda * diag(diag(ld.mat))

            VX <- diag(ld.mat)
            XtY <- n.all * VX * sumstats$BETA
            se <- abs(sumstats$BETA) / -qnorm(sumstats$P/2)
            m1 <- matrix( rep(se^2,k), ncol=k, nrow=k ) *
                matrix( rep(VX,k), ncol=k, nrow=k )
            m <- ifelse( m1<t(m1), m1, t(m1) )
            Sigma <- n.all * m * ld.mat

            m <- XtY * n.prop[1]
            v <- Sigma * n.all * n.prop[1] * (1 - n.prop[1])
            XtY.1 <- mvrnorm( n=n.folds, mu=m, Sigma=v, tol=tol )

            m <- (XtY-t(XtY.1)) * n.prop[2] / (1-n.prop[1])
            v <- Sigma * n.all * n.prop[2] * (1 - n.prop[1]-n.prop[2]) / (1-n.prop[1])
#            XtY.2 <- mvrnorm( n=1, mu=m, Sigma=v )
            XtY.2 <- mvrnorm( n=n.folds, mu=rep(0,k), Sigma=v, tol=tol )
            XtY.2 <- t(t(XtY.2) - m)

            XtY.3 <- t(XtY - t(XtY.1) - t(XtY.2))

            beta.1 <- solve(ld.mat) %*% t(matrix(XtY.1,nrow=n.folds)) / (n.all*n.prop[1])
#            beta.1 <- solve(ld.mat) %*% t(XtY.1) / (n.all*n.prop[1])
            se.1 <- se * sqrt( 1 / n.prop[1] )
            p.1 <- 2*pnorm( abs(beta.1) / se.1, lower.tail=FALSE )
#            sumstats.1 <- data.frame( sumstats$SNP,
#                                     sumstats$ALLELE1, sumstats$ALLELE0,
#                                     beta.1, se.1, p.1, XtY.2, XtY.3 )
            sumstats.1 <- list( sumstats$SNP,
                                     beta.1, se.1, p.1, XtY.2, XtY.3 )
        }else{
#        colnames(sumstats.1) <- c('SNP','ALLELE1','ALLELE0','BETA','SE','P',
#                                  'XtY.2','XtY.3')
            sumstats.1 <- list( NULL,NULL,NULL,NULL,NULL,NULL )
        }
    }else{
        sumstats.1 <- list( NULL,NULL,NULL,NULL,NULL,NULL )
    }
    names(sumstats.1) <- c('SNP','BETA','SE','P','XtY.2','XtY.3')

    return(sumstats.1)
}

read.fit.clump <- function( clump.i, sumstats, ld.ids,
                           X.bed, bim, l=10, S=1, precision=FALSE,
                           by.chr, beta.stem, strand.check ){
    clump.snps <-  unlist(strsplit( clump.i$SP2, ',' ))
    tmp <- strsplit( clump.snps, "\\(1" )
    if ( tmp[[1]][1]!="NONE" ) {
        clump.snps <- c( clump.i$SNP, sapply( tmp, getElement, 1 ) )
    } else {
        clump.snps <- clump.i$SNP
    }
    snps <- intersect( sumstats$SNP, clump.snps )
    chr <- clump.i$CHR
    if( by.chr ){
        X.bed <- X.bed[[chr]]
        bim <- bim[[chr]]
    }

    if( length(snps)>0 ){
        sumstats <- sumstats[match( snps, sumstats$SNP ),]
        ref.stats <- est.ref.stats( snps, ld.ids, X.bed, bim,
                                   sumstats$ALLELE1, sumstats$ALLELE0,
                                   strand.check, n.eff=FALSE )

        ptr.use <- which( !is.na(ref.stats$af) & ref.stats$af!=0 )
        if( length(ptr.use)>0 ){
            ref.stats$af <- ref.stats$af[ptr.use]
            ref.stats$ld <- ref.stats$ld[ptr.use,ptr.use]
            snps <- snps[ptr.use]
            sumstats <- sumstats[ptr.use,]

            n.models <- length(S)*length(l)
            beta.bar <- matrix( ncol=n.models, nrow=length(snps) )
            colnames(beta.bar) <- 1:(n.models)

            if( precision ){
                lambda <- list( length=(n.models) )
            }
            ii <- 0
            for( j in 1:length(l) ){
                for( i in 1:length(S) ){
                    ii <- ii + 1
                    tmp <- ridge.fit( beta=sumstats$BETA,
                                     LD=ref.stats$ld, af=ref.stats$af,
                                     l=l[j], S=S[i], precision=precision )
                    beta.bar[,ii] <- tmp[[1]]
                    colnames(beta.bar)[ii] <- paste('beta.bar',l[j],S[i],sep="_")
#                if( !is.na(tmp[[3]]) ){
#                    names(kl)[ii] <- paste('beta.bar',l[j],S[i],sep="_")
#                }
                    if( precision ){
                        lambda[[ii]] <- tmp[[2]]
                    }
                }
            }
            beta.bar <- data.frame( names(ref.stats$af), clump.i$SNP, clump.i$CHR,
                                   sumstats$ALLELE1, sumstats$ALLELE0,
                                   clump.i$P, ref.stats$af,
                                   sumstats$BETA, beta.bar )
            colnames(beta.bar)[1:8] <- c('snp','clump.id','chr',
                                         'effect.allele','ref.allele',
                                         'p.value','af','beta_hat')
#    rownames(beta.bar) <- names(ref.stats$af)

            ret <- list()
            ret[[1]] <- beta.bar
            if( precision ){
                lambda <- as.data.table( lambda[[1]] )
                fwrite( lambda, paste0(beta.stem,'/lambda/',clump.i$SNP,'.gz'))
            }
            return(ret)
        }
    }
}

read.NonCentralFit.clump <- function( sumstats, ld.ids, X.bed, bim,
                                     beta.prior, lambda.ext, w.prior,
                                     precision=FALSE, by.chr, N.pop2,
                                     ranking, sigma2,
                                     lambda.prior0, S.prior0, strand.check ){
    chr <- beta.prior$chr[1]
    if( by.chr ){
        X.bed <- X.bed[[chr]]
        bim <- bim[[chr]]
    }
    snps <- intersect( bim$V2, intersect( beta.prior$snp, sumstats$SNP ) )
    if( length(snps)>0 ){
        ptr.sumstats <- match( snps, sumstats$SNP )
        ptr.prior <- match( snps, beta.prior$snp )
        sumstats <- sumstats[ptr.sumstats,,drop=FALSE]
        ref.stats <- est.ref.stats( snps, ld.ids, X.bed, bim,
                                   effect.allele=beta.prior$effect.allele[ptr.prior],
                                   ref.allele=beta.prior$ref.allele[ptr.prior],
                                   strand.check, n.eff=(ranking=="pv.minP") )
        ptr.use <- which( 0.005 < ref.stats$af&ref.stats$af < 0.995 )
        if( length(ptr.use) > 0 ){
            if( length(ptr.use) < length(snps) ){
                snps <- snps[ptr.use]
                sumstats <-  sumstats[ptr.use,]
                ptr.prior <- ptr.prior[ptr.use]
                ref.stats$af <- ref.stats$af[ptr.use]
                ref.stats$ld <- ref.stats$ld[ptr.use,ptr.use]
            }
            if( ranking=="pv.ftest" ){
                pv.target <- f.test( sumstats$BETA, ref.stats$ld,
                                    ref.stats$af, N.pop2, sigma2 )
            }
            if( ranking=="thinned.pv.ftest" ){
                pv.target <- f.test.diag( sumstats$BETA, ref.stats$ld,
                                         ref.stats$af, N.pop2, sigma2 )
            }
            if( ranking=="pv.minP" ){
                pv.target <- min(sumstats$P) * ref.stats$n.eff
            }
            infile <- paste0( lambda.ext, beta.prior$clump.id[1],'.gz' )
            lambda.prior <- as.matrix(fread(infile))

            if( length(ptr.prior) == nrow(beta.prior) ){
                lambda.prior <- lambda.prior[ptr.prior,ptr.prior,drop=FALSE]
                beta.prior <- beta.prior[ptr.prior,]
            }
            if( length(ptr.prior) < nrow(beta.prior) ){
                beta <- as.matrix(beta.prior$beta.bar)
                xTy <- ( lambda.prior %*% beta )[ptr.prior,,drop=FALSE]
                lambda.prior <- lambda.prior[ptr.prior,ptr.prior,drop=FALSE]
                beta <- solve(lambda.prior) %*% xTy
                beta.prior <- beta.prior[ptr.prior,]
                beta.prior$beta.bar <- beta
            }

#        beta.prior <- beta.prior[ptr.prior,,drop=FALSE]
#        sigma.prior <- solve(lambda.prior)[ptr.prior,ptr.prior,drop=FALSE]
#        lambda.prior <- solve(sigma.prior)

            af.pop1 <- beta.prior$af
            var.pop1 <- 2*af.pop1*(1-af.pop1)
            lambda.prior0 <- lambda.prior0 * var.pop1^S.prior0

            beta.bar <- matrix( ncol=length(w.prior), nrow=length(snps) )
            colnames(beta.bar) <- 1:(length(w.prior))
            kl <- vector(length=length(w.prior))

            swtch <- allele.check( beta.prior$effect.allele, beta.prior$ref.allele,
                                  sumstats$ALLELE1, sumstats$ALLELE0, strand.check )
#            swtch <- ifelse( beta.prior$effect.allele==sumstats$ALLELE1 &
#                             beta.prior$ref.allele==sumstats$ALLELE0, 1, 0 )
#            swtch <- ifelse( beta.prior$effect.allele==sumstats$ALLELE0 &
#                             beta.prior$ref.allele==sumstats$ALLELE1, -1, swtch )
#            ptr.miss <- which( swtch==0 )
#            if( strand.check & length(ptr.miss)>0 ){
#                coded.allele1 <- alt.strand( sumstats$ALLELE1 )
#                other.allele1 <- alt.strand( sumstats$ALLELE0 )
#                swtch[ptr.miss] <- ifelse( coded.allele1[ptr.miss]==
#                                           beta.prior$effect.allele[ptr.miss] &
#                                           other.allele1[ptr.miss]==
#                                           beta.prior$ref.allele[ptr.miss],
#                                          1, 0 )
#                swtch[ptr.miss] <- ifelse( coded.allele1[ptr.miss]==
#                                           beta.prior$ref.allele[ptr.miss] &
#                                           other.allele1[ptr.miss]==
#                                           beta.prior$effect.allele[ptr.miss],
#                                          -1, swtch[ptr.miss] )
#            }

            # Map sumstats to prior
            if( precision ){
                lambda <- list( length=length(w.prior) )
            }
            for( i in 1:length(w.prior) ){
                tmp <- noncentral.ridge.fit( beta.data=sumstats$BETA*swtch,
                                            LD=ref.stats$ld, af=ref.stats$af,
                                            beta.tilde1=beta.prior$beta.bar,
                                            lambda0=lambda.prior0,
                                            lambda1=lambda.prior,
                                            w.prior=w.prior[i],
                                            n=N.pop2, precision=precision,
                                            ranking=ranking, sigma2=sigma2 )
                beta.bar[,i] <- tmp[[1]]
                colnames(beta.bar)[i] <- paste('beta.bar',w.prior[i],sep="_")
                if( ranking=="f.stat" | ranking=="thinned.f.stat" ){
                    kl[i] <- tmp[[3]]
                    names(kl)[i] <- paste('beta.bar',w.prior[i],sep="_")
                }
                if( precision ){
                    lambda[[i]] <- tmp[[2]]
                }
            }

            if( ranking=="pv.ftest" | ranking=="thinned.pv.ftest" |
                ranking=="pv.minP" ){
                beta.bar <- data.frame( beta.prior[,1:6],
                                       sumstats$BETA, pv.target, beta.prior$beta.bar,
                                       beta.bar )
                colnames(beta.bar)[1:9] <- c('snp','clump.id','chr',
                                             'effect.allele','ref.allele','p.value',
                                             'beta_hat','pv.target','beta_prior')
            }else{
                beta.bar <- data.frame( beta.prior[,1:6],
                                       sumstats$BETA, beta.prior$beta.bar,
                                       beta.bar )
                colnames(beta.bar)[1:8] <- c('snp','clump.id','chr',
                                             'effect.allele','ref.allele','p.value',
                                             'beta_hat','beta_prior')
            }
            rownames(beta.bar) <- snps

            ret <- list(length=3)
            ret[[1]] <- beta.bar
            ret[[2]] <- kl
            if( precision ){
                names(lambda) <- colnames(beta.bar)[-(1:9)]
                fwrite( lambda, paste0(beta.stem,'/lambda/',clump.i$SNP,'.gz') )
#            ret[[2]] <- lambda
            }
            return(ret)
        }
    }
}

get.pred.clump <- function( beta.bar, ptr.beta.use, clump.id, X.bed, bim,
                           ids.use, strand.check ){
    snps <- intersect( beta.bar$snp, bim$V2 )
    if(length(snps) < length(beta.bar$snp)){
        print("Warning!: intersection is lower than beta.bar!")
        print(paste0(length(snps)," > ", length(beta.bar$snp)))
        beta.bar <- beta.bar[beta.bar$snp %in% snps,]
    }
    X <- X.bed[ ids.use, beta.bar$snp, drop=FALSE ]

    ptr <- match( beta.bar$snp, colnames(X) )
    X <- X[,ptr,drop=FALSE]
    ptr.bed <- match( beta.bar$snp, colnames(X.bed) )

    coded.allele <- bim$V5[ptr.bed]
    other.allele <- bim$V6[ptr.bed]

    swtch <- allele.check( coded.allele, other.allele,
                          beta.bar$effect.allele, beta.bar$ref.allele, strand.check )
#    swtch <- ifelse( coded.allele==beta.bar$effect.allele &
#                     other.allele==beta.bar$ref.allele, 1, 0 )
#    swtch <- ifelse( coded.allele==beta.bar$ref.allele &
#                     other.allele==beta.bar$effect.allele, -1, swtch )
#    ptr.miss <- which( swtch==0 )
#    if( strand.check & length(ptr.miss)>0 ){
#        coded.allele1 <- alt.strand( coded.allele )
#        other.allele1 <- alt.strand( other.allele )
#        swtch[ptr.miss] <- ifelse( coded.allele1[ptr.miss] ==
#                                   beta.bar$effect.allele[ptr.miss] &
#                                   other.allele1[ptr.miss] ==
#                                   beta.bar$ref.allele[ptr.miss],
#                                  1, 0 )
#        swtch[ptr.miss] <- ifelse( coded.allele1[ptr.miss]==
#                                   beta.bar$ref.allele[ptr.miss] &
#                                   other.allele1[ptr.miss]==
#                                   beta.bar$effect.allele[ptr.miss],
#                                  -1, swtch[ptr.miss] )
#    }

#    Map target genotypes to estimated SNP effects
    ptr.use <- which( swtch!=0 )
    if( length(ptr.use)>0 ){
        X <- X[,ptr.use,drop=FALSE]
        beta.bar <- beta.bar[ptr.use,,drop=FALSE]

        ptr <- which(swtch == -1)
        if( length(ptr)>0 ){
            X[,ptr] <- 2 - X[,ptr]
        }

        m <- apply( X, 2, mean, na.rm=TRUE )
        for( ii in 1:ncol(X) ){
            X[,ii] <- ifelse( is.na(X[,ii]), m[ii], X[,ii] )
        }
        pred <- as.matrix(X) %*% as.matrix(beta.bar[,ptr.beta.use])
        ret <- pred
    }else{
        ret <- matrix(data=0,ncol=length(ptr.beta.use),nrow=nrow(X))
    }
    rownames(ret) <- rownames(X)

    return(ret)
}

get.pred.genome <- function( beta.bar, p.thresh, X.bed, bim,
                            ids.use, b.size=100, n.cores=20,
                            kl.metric, kl.thresh, ranking, strand.check ){
    n.clumps <- length(beta.bar)
    batches <- ceiling(n.clumps/b.size)
    chr <- as.numeric(sapply( sapply(beta.bar,getElement,'chr'), getElement, 1 ))
    if( ranking=="pv" ){
        pval <- as.numeric(sapply( sapply(beta.bar,getElement,'p.value'),
                                  getElement, 1 ))
    }
    if( ranking=="pv.minP" | ranking=="pv.ftest" | ranking=="thinned.pv.ftest" ){
        pval <- as.numeric(sapply( sapply(beta.bar,getElement,'pv.target'),
                                  getElement, 1 ))
    }

    if( ranking=="f.stat" | ranking=="thinned.f.stat" ){
        pred.genome <- list(length=nrow(kl.thresh))
        ptr.beta.use <- match( colnames(kl.metric), colnames(beta.bar[[1]]) )
        nc <- length(ptr.beta.use)
        for( j in 1:nrow(kl.thresh) ){
            pred.genome[[j]] <- matrix( ncol=nc, nrow=length(ids.use), data=0 )
            rownames(pred.genome[[j]]) <- ids.use
            colnames(pred.genome[[j]]) <- colnames(kl.metric)
        }
    }
    if( ranking=="pv" | ranking=="pv.minP" |
        ranking=="pv.ftest" | ranking=="thinned.pv.ftest" ){
        pred.genome <- list(length=nrow(p.thresh))
        ptr.beta.use <- grep( 'beta.bar', colnames(beta.bar[[1]]) )
        nc <- length(ptr.beta.use)
        for( j in 1:nrow(p.thresh) ){
            pred.genome[[j]] <- matrix( ncol=nc, nrow=length(ids.use), data=0 )
        }
    }

    for( k in 1:batches ){
        mx <- ifelse( k*b.size>n.clumps, n.clumps, k*b.size )
        i.clump.use <- (1+(k-1)*b.size):mx
        pred <- mclapply( i.clump.use,
                         function(i){
                             get.pred.clump( beta.bar=beta.bar[[i]],
                                            ptr.beta.use=ptr.beta.use,
                                            clump.id=names(beta.bar)[i],
                                            X.bed=X.bed, bim=bim,
                                            ids.use=ids.use,
                                            strand.check=strand.check )},
                         mc.cores=n.cores )
#        for( i in i.clump.use ){
#            pred <- get.pred.clump( beta.bar=beta.bar[[i]], ptr.beta.use=ptr.beta.use, clump.id=names(beta.bar)[i], X.bed=X.bed, bim=bim, by.chr=by.chr, ids.use=ids.use, strand.check=strand.check )
#            print(c(i,range(pred)))
#        }
        if( ranking=='f.stat' | ranking=='thinned.f.stat' ){
            for( j in 1:nrow(kl.thresh) ){
                for( kk in 1:nc ){
                    ptr.use <- which( kl.metric[i.clump.use,kk] >= kl.thresh[j,kk] )
                    for( i in ptr.use ){
                        pred.genome[[j]][,kk] <- pred.genome[[j]][,kk] + pred[[i]][,kk]
                    }
                }
            }
        }else{
            for( j in 1:nrow(p.thresh) ){
                ptr.use <- which( p.thresh[j,1] >= pval[i.clump.use]&
                                  pval[i.clump.use] > p.thresh[j,2] )
                for( i in ptr.use ){
                    pred.genome[[j]] <- pred.genome[[j]] + pred[[i]]
                }
            }
        }
    }
    return(pred.genome)
}

thin.big.loci <- function( clump, thinned.snplist, n.max.locus, sumstats2=NULL ){
    for( i in 1:nrow(clump) ){
        clump.snps <-  unlist(strsplit( clump$SP2[i], ',' ))
        tmp <- strsplit( clump.snps, "\\(1" )
        clump.snps <- sapply( tmp, getElement, 1 )
        if( length(clump.snps)>n.max.locus ){
            clump.snps.thinned <- intersect( clump.snps, thinned.snplist )
            if( !is.null(sumstats2) ){
                snps.pop2 <- intersect( clump.snps, sumstats2$SNP )
                ptr.pop2 <- match( snps.pop2, sumstats2$SNP )
                s <- order( sumstats2$P[ptr.pop2], decreasing=FALSE )
                clump.snps.thinned <-  c( sumstats2$SNP[ptr.pop2[s[1]]], clump.snps.thinned )
                clump.snps.thinned <- unique(clump.snps.thinned)
            }
            clump$SP2[i] <- paste0( clump.snps.thinned, collapse='(1),' )
        }
    }
    return(clump)
}

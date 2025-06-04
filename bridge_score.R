library(data.table)

dir <- 'out/preds/'
model <- c('best','weighted')
pop <- c('AFR','EAS')
pred <- list()
pop <- 'EAS'
pheno <- fread(paste0('data/pop_',pop,'/phenotypes/',pop,'_valid.dat'))
for( i in 1:2 ){
    tmp <- fread(paste0(dir,pop,'_',model[i],'_chr1.profile'))
    pred[[i]] <- data.frame( tmp$IID, tmp$SCORE * tmp$CNT )
    for( chr in 2:22 ){
        tmp <- fread(paste0(dir,pop,'_',model[i],'_chr',chr,'.profile'))
        pred[[i]][,2] <- pred[[i]][,2] + tmp$SCORE * tmp$CNT
    }
    print(paste(model[i],cor( pheno$y, pred[[i]][,2] )^2))
}

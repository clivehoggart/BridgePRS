library(data.table)

dir <- 'out/preds/'
model <- c('best','weighted')
pop <- c('AFR','EAS')
pred <- list()
pop <- 'EAS'
pheno <- fread(paste0('data/pop_',pop,'/phenotypes/',pop,'_valid.dat'))
for( i in 1:2 ){
    tmp <- fread(paste0(dir,pop,'_',model[i],'_chr1.sscore'))
    pred[[i]] <- data.frame( tmp$IID, tmp$SCORE1_AVG*tmp$ALLELE_CT )
    for( chr in 2:22 ){
        tmp <- fread(paste0(dir,pop,'_',model[i],'_chr',chr,'.sscore'))
        pred[[i]][,2] <- pred[[i]][,2] + tmp$SCORE1_AVG * tmp$ALLELE_CT
    }
    print(paste(model[i],cor( pheno$y, pred[[i]][,2] )^2))
}

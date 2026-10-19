#BSUB -L /bin/sh
#BSUB -n 20
#BSUB -J Bridge[1-3,5-6,8,13-18,21-24,28-33,37-84,91-99,103-114,121-132,136-141,145-147,151-153,157-162,166-175,177-190]
#BSUB -R "span[hosts=1]"
#BSUB -q premium
#BSUB -W 48:00
#BSUB -P acc_psychgen
#BSUB -o Bridge.o%J.%I
#BSUB -eo Bridge.e%J.%I
#BSUB -M 20000

module load R/4.5.1
module load plink

rpt=10
idx=${LSB_JOBINDEX}
phenos=(50 21001 30080 30710 30140 30040 30630 30240 30610 30780 30130 30070 30670 30870 30220 30860 30530 30770 30150) 
i=$(( (idx - 1) / rpt ))
iter=$(( (idx - 1) % rpt + 1 ))
pheno=${phenos[$i]} 

echo "Trait: $pheno"
echo "Run:   $iter"

home=/hpc/users/hoggac01
data=$home/BridgePRS/data
src=$home/BridgePRS/src

pop1=EUR
pop2=AFR
pop22=afr
fst=0.15
outdir=/sc/arion/projects/psychgen/projects/prs/cross_population_prs_development/quick_ridge/results/sumstats/imputed
path=/sc/arion/projects/psychgen/ukb/usr/clive/ukb

#    --by_chr_sumstats _$pop2.gz \
#    --by_chr_sumstats .Phenotype.glm.linear.gz \

${src}/Bash/BridgePRS_sumstat_parallel.sh $src/Rscripts \
    --strand_check 0 \
    --outdir $outdir/$pheno \
    --blockdir $outdir/blocks \
    --by_chr 1 \
    --pop1 $pop1 \
    --pop2 $pop2 \
    --fst $fst \
    --N_pop1 250000 \
    --N_pop2 4800 \
    --pop1_qc_snplist $path/imputed/chr.EUR.qced.snps.snplist \
    --pop2_qc_snplist $path/imputed/chr.$pop2.qced.snps.snplist \
    --pop1_sumstats   $path/impute_gwas/$pheno/$pop1.chr \
    --pop2_sumstats   $path/impute_gwas/$pheno/$pop2.chr \
    --by_chr_sumstats .Phenotype.glm.linear.gz \
    --pop1_ld_bfile   /sc/arion/projects/paul_oreilly/data/Biobanks/UKB/imputed/qced_data/chr \
    --pop1_ld_ids     $path/ld_ids.txt \
    --pop2_ld_bfile   /sc/arion/projects/paul_oreilly/data/Biobanks/UKB/imputed/qced_data/chr \
    --pop2_ld_ids     $path/${pop22}_train_ids.txt \
    --sumstats_snpID ID \
    --sumstats_p P \
    --sumstats_beta BETA \
    --sumstats_allele1 A1 \
    --sumstats_allele0 AX \
    --pheno_name y \
    --n_cores 20 \
    --prop_train 0.6 \
    --prop_test 0.3 \
    --do_sumstat_pop1 1 \
    --do_clump_pop1 1 \
    --do_est_beta_pop1 1 \
    --do_sumstat_ensembl_pop1 1 \
    --do_est_beta_pop1_precision 1 \
    --do_sumstat_pop2 1 \
    --do_est_beta_InformPrior 1 \
    --do_clump_pop2 1 \
    --do_est_beta_pop2 1 \
    --do_sumstat_ensembl_pop2 1 \
    --clean 0 \
    --iter $iter


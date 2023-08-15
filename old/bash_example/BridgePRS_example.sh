#!/bin/bash

#DATAPATH=~/BridgePRS1/data
AFR_DATA=../test_data/pop1_africa
EUR_DATA=../test_data/pop2_europe
BFILE_PREFIX=../test_data/plink_files/bfile3  


./BridgePRS.sh   \
    --strand_check 1 \
    --outdir out \
    --by_chr 0 \
    --by_chr_sumstats .Phenotype.glm.linear.gz \
    --pop1 EUR \
    --pop2 AFR \
    --fst 0.15 \
    --pop1_sumstats   ${EUR_DATA}/sumstats/EUR.chr \
    --pop1_qc_snplist ${EUR_DATA}/snps.txt \
    --pop1_ld_ids     ${EUR_DATA}/ids.txt \
    --pop1_test_data  ${EUR_DATA}/phenos.dat \
    --pop2_sumstats   ${AFR_DATA}/sumstats/AFR.chr \
    --pop2_qc_snplist ${AFR_DATA}/snps.txt \
    --pop2_ld_ids     ${AFR_DATA}/ids.txt \
    --pop2_test_data  ${AFR_DATA}/phenos.test.dat \
    --pop2_valid_data ${AFR_DATA}/phenos.valid.dat \
    --pop1_bfile      $BFILE_PREFIX \
    --pop1_ld_bfile   $BFILE_PREFIX \
    --pop2_bfile      $BFILE_PREFIX \
    --pop2_ld_bfile   $BFILE_PREFIX \
    --sumstats_snpID ID \
    --sumstats_p P \
    --sumstats_beta BETA \
    --sumstats_allele1 A1 \
    --sumstats_allele0 AX \
    --sumstats_n OBS_CT \
    --sumstats_se SE \
    --sumstats_frq A1_FREQ \
    --pheno_name y \
    --cov_names PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --n_cores 6 \
    --do_clump_pop1 1 \
    --do_est_beta_pop1 1 \
    --do_predict_pop1 1 \
    --do_est_beta_pop1_precision 1 \
    --do_est_beta_InformPrior 1 \
    --do_predict_pop2_stage2 1 \
    --do_clump_pop2 1 \
    --do_est_beta_pop2 1 \
    --do_predict_pop2 1 \
    --do_combine 1 \

#!/bin/sh

data="$PWD"/data 
src="$PWD"/src 

pop=(AFR EAS)
fst=(0.15 0.11)
i=1
${src}/Bash/BridgePRS_sumstat.sh $src/Rscripts \
    --strand_check 1 \
    --outdir out \
    --by_chr 1 \
    --by_chr_sumstats .glm.linear.gz \
    --pop1 EUR \
    --pop2 ${pop[$i]} \
    --fst ${fst[$i]} \
    --N_pop1 80000 \
    --N_pop2 20000 \
    --pop1_qc_snplist ${data}/qc_snplist.txt \
    --pop1_sumstats   ${data}/pop_EUR/sumstats/EUR.chr \
    --pop2_qc_snplist ${data}/qc_snplist.txt \
    --pop2_sumstats   ${data}/pop_${pop[$i]}/sumstats/${pop[$i]}.chr \
    --pop1_ld_bfile   ${data}/1000G_sample/chr \
    --pop1_ld_ids     ${data}/1000G_sample/EUR_IDs.txt \
    --pop2_ld_bfile   ${data}/1000G_sample/chr \
    --pop2_ld_ids     ${data}/1000G_sample/${pop[$i]}_IDs.txt \
    --sumstats_snpID ID \
    --sumstats_p P \
    --sumstats_beta BETA \
    --sumstats_allele1 ALT \
    --sumstats_allele0 REF \
    --pheno_name y \
    --n_cores 20 \
    --prop_train 0.6 \
    --prop_test 0.3 \
    --do_block_pop1 1 \
    --do_block_pop2 1 \
    --do_sumstat_pop1 1 \
    --do_sumstat_pop2 1 \
    --do_clump_pop1 1 \
    --do_est_beta_pop1 1 \
    --do_sumstat_ensembl_pop1 1 \
    --do_est_beta_pop1_precision 1 \
    --do_est_beta_InformPrior 1 \
    --do_clump_pop2 1 \
    --do_est_beta_pop2 1 \
    --do_sumstat_ensembl_pop2 1 \
    --do_pool 1 \
    --n_folds 2

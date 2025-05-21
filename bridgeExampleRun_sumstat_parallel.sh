#BSUB -L /bin/sh
#BSUB -n 8
#BSUB -J Bridge[1-10]
#BSUB -R "span[hosts=1]"
#BSUB -q premium
#BSUB -W 1:00
#BSUB -P acc_psychgen
#BSUB -o Bridge.o%J.%I
#BSUB -eo Bridge.e%J.%I
#BSUB -M 20000

module load R/4.1.0
module load plink

iter=${LSB_JOBINDEX}

data="$PWD"/data 
src="$PWD"/src 

pop=(AFR EAS)
fst=(0.15 0.11)
i=1
${src}/Bash/BridgePRS_sumstat_parallel.sh $src/Rscripts \
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
    --n_cores 8 \
    --prop_train 0.6 \
    --prop_test 0.3 \
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
    --fold $iter

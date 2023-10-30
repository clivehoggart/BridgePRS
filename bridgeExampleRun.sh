data="$PWD"/data 
src="$PWD"/src 

${src}/Bash/BridgePRS.sh $src/Rscripts\
    --strand_check 1 \
    --outdir out_binary \
    --by_chr 1 \
    --by_chr_sumstats .glm.linear.gz \
    --pop1 EUR \
    --pop2 AFR \
    --fst 0.15 \
    --pop1_sumstats   ${data}/sumstats/eur.chr \
    --pop1_qc_snplist ${data}/sumstats/qc_snplist3.txt \
    --pop1_bfile      ${data}/bfiles/chr \
    --pop1_test_data  ${data}/eur_test.dat \
    --pop1_ld_bfile   ${data}/bfiles/chr \
    --pop1_ld_ids     ${data}/eur_ids.txt \
    --pop2_sumstats   ${data}/sumstats/afr.chr \
    --pop2_qc_snplist ${data}/sumstats/qc_snplist3.txt \
    --pop2_bfile      ${data}/bfiles/chr \
    --pop2_test_data  ${data}/afr_test.dat \
    --pop2_valid_data ${data}/afr_valid.dat \
    --pop2_ld_bfile  ${data}/bfiles/chr \
    --pop2_ld_ids    ${data}/afr_ids.txt \
    --sumstats_snpID ID \
    --sumstats_p P \
    --sumstats_beta BETA \
    --sumstats_allele1 ALT \
    --sumstats_allele0 REF \
    --sumstats_n OBS_CT \
    --sumstats_se SE \
    --sumstats_frq A1_FREQ \
    --pheno_name y \
    --n_cores 20 \
    --do_clump_pop1 1 \
    --do_est_beta_pop1 1 \
    --do_predict_pop1 1 \
    --do_est_beta_pop1_precision 1 \
    --do_est_beta_InformPrior 1 \
    --do_predict_pop2_stage2 1 \
    --do_clump_pop2 1 \
    --do_est_beta_pop2 1 \
    --do_predict_pop2 1 \
    --do_combine 1 


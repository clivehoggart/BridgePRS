data="$PWD"/data 
src="$PWD"/src 

pop=(AFR EAS)
fst=(0.15 0.11)
i=0
${src}/Bash/BridgePRS.sh $src/Rscripts\
    --strand_check 1 \
    --outdir out \
    --by_chr 1 \
    --by_chr_sumstats .glm.linear.gz \
    --pop1 EUR \
    --pop2 AFR \
    --fst ${fst[$i]} \
    --pop1_qc_snplist ${data}/qc_snplist.txt \
    --pop1_sumstats   ${data}/pop_EUR/sumstats/EUR.chr \
    --pop1_bfile      ${data}/pop_EUR/genotypes/chr \
    --pop1_test_data  ${data}/pop_EUR/phenotypes/EUR_test.dat \
    --pop2_qc_snplist ${data}/qc_snplist.txt \
    --pop2_sumstats   ${data}/pop_${pop[$i]}/sumstats/${pop[$i]}.chr \
    --pop2_bfile      ${data}/pop_${pop[$i]}/genotypes/chr \
    --pop2_test_data  ${data}/pop_${pop[$i]}/phenotypes/${pop[$i]}_test.dat \
    --pop2_valid_data ${data}/pop_${pop[$i]}/phenotypes/${pop[$i]}_valid.dat \
    --pop1_ld_bfile   ${data}/1000G_sample/chr \
    --pop1_ld_ids     ${data}/1000G_sample/EUR_IDs.txt \
    --pop2_ld_bfile   ${data}/1000G_sample/chr \
    --pop2_ld_ids     ${data}/1000G_sample/${pop[$i]}_IDs.txt \
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

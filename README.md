# BridgePRS
Estimation of trans-ancestry PRS 
BridgePRS is an R and bash based package that integrates GWAS summary
statistics from two populations. It was designed to improve prediction
in a population for which GWASs are relatively under-powered
(population 2) but there exist powerful GWAS data in another
population (population 1). In addition to GWAS summary statistics,
BridgePRS also requires genotype and phenotype data from the two
populations for parameter optimisation (test data) and estimation of
LD. The LD and test data files can be different or identical. Measures
of model fit and SNP weights for the best model for population 2 are
returned. In addition, validation data from population 2, can be
supplied to make out-of-sample prediction and measures of model fit.

# Installation
Clone this repository using the following git command:

git clone https://github.com/clivehoggart/BridgePRS.git

Alternatively, download the source files from the github website
(https://github.com/clivehoggart/BridgePRS.git)

Required R packages: parallel, glmnet, boot, data.table, doMC,
BEDMatrix v2.0.3, optparse and MASS.
Assumes Plink v1.9 is in your path

# Getting started
Example script can be run from the installation directory by typing.
./Bridge_example.sh
Script utilises 8 processors for paralel computation, this can be changed
from the --n_cores argument, see below.

BridgePRS is implemented by clumping GWAS summary statistics from
populations 1 and 2 in plink using the following clumping arguments

--clump-p1 1e-1 --clump-p2 1e-1 --clump-kb 1000 --clump-r2 0.01

and subsequently running a series of R scripts to estimate PRS and
tune parameters by fitting the PRS to test data. The individual steps
can be run separately, in particular code will run faster if user
parallelises clumping step by running by chromosomes.

# To run BridgePRS:
BridgePRS/bin/ridgePRS_bychr.sh \
	--outdir. Output directory
	
	--pop1. Label for population 1 result files (default pop1)
	
	--pop2. Label for population 2 result files (default pop2)

	--pop1_sumstats. Population 1 sumstats file, can be gzipped (required)

	--pop2_sumstats. Population 2 sumstats file, can be gzipped (required)

	--by_chr_sumstats. If set assumes sumstats are split by
          chromosome and are labelled pop1_sumstatsXby_chr_sumstats,
          where X indexes the chromosomes

	--pop1_qc_snplist. List of variants to use in clumping of
          population 2 sumstats (required)

	--pop2_qc_snplist. As above for population 2 (required).

	--pop1_bfile. Population 1 plink binary file to calculate PRS
          in test and validation sets (required)

	--pop1_test_data. Population 1 phenotype and covariate test
          file. Requires a column labelled IID of IDs to use
          (required)

	--pop1_valid_data. As above but for validation data (optional)

	--pop1_ld_bfile. Population 1 plink binary file for LD
	  estimation (required)

	--pop1_ld_ids. Population 1 IDs of individuals to use for LD
	  estimation (required)

	--pop2_bfile. As above for population 2 (required)

	--pop2_test_data. As above for population 2 (required)

	--pop2_valid_data. As above for population 2 (optional)

	--pop2_ld_bfile. As above for  population 2 (required)

	--pop2_ld_ids. As above for population 2 (required)

	--cov_names. Column names of covariates to control for in
	  optimising PRS and calculating R2 (optional)

	--pheno_name. Column name of phenotype used in test and
	  validation files (required)

    	--sumstats.snpID. Sumstat column label of SNP ID (required)

    	--sumstats.betaID. Sumstat column label of SNP effect (required)

    	--sumstats.allele1ID. Sumstat column label SNP effect allele (required)

   	--sumstats.allele0ID. Sumstat column label of SNP ref allele (required)

	--n_cores. Nummber of processors to use by R scripts

	--do_clump_pop1. Logical to run plink to clump population 1
          summary statistis (default 1)
    
	--do_est_beta_pop1. Logical to run R script to estimate PRS using
          population 1 summary statistics (default 1)

    	--do_predict_pop1. Logical to run R script to make predictions
          in population 1 in test and if supplied validation data
          (default 1)
    
	--do_est_beta_pop1_precision. Logical to run R script to
          output precision matrices of best fitting model in
          population 1 (default 1)

	--do_est_beta_InformPrior. Logical to run R script to estimate
          PRS in population 2 given prior from population 1 (default 1)

	--do_predict_pop2_stage2. Logical to run R script to make
          prediction in population 2 in test and if supplied
          validation data using popultion 1 informed PRS (default 1)

	--do_clump_pop2. Logical to run plink to clump population 2
          summary statistis (default 1)

	--do_est_beta_pop2. Logical to run R script to estimate PRS
          using population 2 summary statistics (default 1)

	--do_predict_pop2. Logical to run R script to to make
          predictions in population 2 in test and if supplied
          validation data using PRS informed only by population 2
          summary statistics (default 1)

	--do_combine. Logical to run R script to combine PRS in
          population 2. Outputs R2 in test and validation data and SNP
          weights of best model (default 1)

# Output files

pop_stageX_all_preds_test.dat
	-- PRS predictions for test individuals for prior parameters
pop_stageX_all_preds_valid.dat    
	-- PRS predictions for validation individuals for prior parameters
pop_stageX_best_model_params.dat
	-- best fitting model parameters in test set
pop_stageX_best_pred_test.dat     
	-- PRS predictions for test individuals for best prior parameters
pop_stageX_best_pred_valid.dat    
	-- PRS predictions for valid individuals for best prior parameters

X=1 -- PRS estimated using sumstats of target population
X=2 -- PRS estimated using using sumstats from both populations

AFR_weighted_combined_preds.dat
	--  PRS predictions in validation samples of best weighted PRS
	    across all estimated PRS 
AFR_weighted_combined_snp_weights.dat
	--  PRS SNP weights of best weighted PRS across all estimated
	    PRS with columns SNP ID, effect allele, ref allele and effect.
AFR_weighted_combined_var_explained.txt
	-- Variance explained of best weighted PRS in validation data

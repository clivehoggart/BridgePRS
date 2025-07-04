#!/bin/bash
# update key
# Set some default values:

safe_rm() {
    rm "$@" 2>/dev/null || {
        err=$?
        if [ $err -ne 0 ] && [ $err -ne 1 ]; then
            echo "ERROR: rm failed on $*"
            exit $err
        fi
    }
}

set -e
trap 'echo "ERROR on line $LINENO: $BASH_COMMAND"; exit 1' ERR

RSCRIPTS=~/BridgePRS/src/Rscripts
FPATH=$RSCRIPTS"/functions.R"

pheno_name="y"
ext=""
n_cores=1

do_block_pop1=0
do_block_pop2=0
do_sumstat_pop1=0
do_sumstat_pop2=0
do_clump_pop1=0
do_clump_pop2=0
do_est_beta_pop1=0
do_sumstat_ensembl_pop1=0
do_est_beta_pop1_precision=0
do_est_beta_InformPrior=0
do_est_beta_pop2=0
do_sumstat_ensembl_pop1=0
do_pool=0

by_chr=0
by_chr_ld=0
by_chr_target=0
by_chr_sumstats=0
pop1="pop1"
pop2="pop2"
thinned_snplist=0
n_max_locus=0
ids_col=TRUE
ranking=f.stat
pop1_valid_data=0
pop2_valid_data=0
strand_check=0
pop1_ld_bfile=0
pop2_ld_bfile=0
cov_names="000"

n_folds=5
prop_train=0.7
prop_test=0.15
blockdir=0

clean=1

#RSCRIPTS=~/BridgePRS/src/Rscripts
#FPATH=$RSCRIPTS"/functions.R"


#usage()
#{
#  echo "Usage: ridgePRS [ -n | --n_cores n_cores ] [ -b | --bfile bfile ]
#                        [ -o | --outdir outdir ] 
#                        [ -v | --eur_ld_ids eur_ld_ids ] filename(s)"
#  exit 2
#}

PARSED_ARGUMENTS=$(getopt -a -n ridgePRS -o b:o:n:c:d:e:f:g:h:i:j:k:l:m:p:q:r:s:t:u:v:w:x:y:z:1:2:3:4:5:6:7: --long outdir:,blockdir:,n_cores:,pop1_ld_ids:,pop2_ld_ids:,pop1_ld_bfile:,pop2_ld_bfile:,pop1_sumstats:,pop2_sumstats:,pop1_valid_data:,pop2_valid_data:,pop1_test_data:,pop2_test_data:,pop1_qc_snplist:,pop2_qc_snplist:,do_clump_pop1:,do_est_beta_pop1:,,do_est_beta_pop1_precision:,do_est_beta_InformPrior:,do_clump_pop2:,do_est_beta_pop2:,do_est_beta_pop2:,do_block_pop1:,do_sumstat_pop1:,do_block_pop2:,do_sumstat_pop2:,do_sumstat_ensembl_pop1:,do_sumstat_ensembl_pop2:,by_chr:,by_chr_target:,by_chr_ld:,cov_names:,pheno_name:,by_chr_sumstats:,pop2:,thinned_snplist:,n_max_locus:,ranking:,pop1:,ids_col:,sumstats_snpID:,sumstats_beta:,sumstats_allele1:,sumstats_allele0:,sumstats_p:,sumstats_n:,sumstats_se:,sumstats_frq:,strand_check:,fst:,binary:,N_pop1:,N_pop2:,do_pool:,n_folds:,prop_train:,prop_test:,clean: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

#echo "PARSED_ARGUMENTS is $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    --outdir) outdir="$2" ; shift 2 ;;
    --blockdir) blockdir="$2" ; shift 2 ;;
    --n_cores) n_cores="$2" ; shift 2 ;;
    --pop1_ld_ids) pop1_ld_ids="$2" ; shift 2 ;;
    --pop2_ld_ids) pop2_ld_ids="$2" ; shift 2 ;;
    --pop1_ld_bfile) pop1_ld_bfile="$2" ; shift 2 ;;
    --pop2_ld_bfile) pop2_ld_bfile="$2" ; shift 2 ;;
    --pop1_sumstats)  pop1_sumstats="$2" ; shift 2  ;;
    --pop2_sumstats)  pop2_sumstats="$2"   ; shift 2  ;;
    --pop1_valid_data) pop1_valid_data="$2" ; shift 2 ;;
    --pop2_valid_data) pop2_valid_data="$2" ; shift 2 ;;
    --pop1_test_data) pop1_test_data="$2" ; shift 2 ;;
    --pop2_test_data) pop2_test_data="$2" ; shift 2 ;;
    --pop1_qc_snplist) pop1_qc_snplist="$2" ; shift 2 ;;
    --pop2_qc_snplist) pop2_qc_snplist="$2" ; shift 2 ;;
    --do_block_pop1) do_block_pop1="$2" ; shift 2 ;;
    --do_sumstat_pop1) do_sumstat_pop1="$2" ; shift 2 ;;
    --do_block_pop2) do_block_pop2="$2" ; shift 2 ;;
    --do_sumstat_pop2) do_sumstat_pop2="$2" ; shift 2 ;;
    --do_clump_pop1) do_clump_pop1="$2" ; shift 2 ;;
    --do_est_beta_pop1) do_est_beta_pop1="$2" ; shift 2 ;;
    --do_sumstat_ensembl_pop1) do_sumstat_ensembl_pop1="$2" ; shift 2 ;;
    --do_est_beta_pop1_precision) do_est_beta_pop1_precision="$2" ; shift 2 ;;
    --do_est_beta_InformPrior) do_est_beta_InformPrior="$2" ; shift 2 ;;
    --do_sumstat_ensembl_pop2) do_sumstat_ensembl_pop2="$2" ; shift 2 ;;
    --do_clump_pop2) do_clump_pop2="$2" ; shift 2 ;;
    --do_est_beta_pop2) do_est_beta_pop2="$2" ; shift 2 ;;
    --do_pool) do_pool="$2" ; shift 2 ;;
    --by_chr) by_chr=$2 ; shift 2 ;;
    --by_chr_target) by_chr_target=$2 ; shift 2 ;;
    --by_chr_ld) by_chr_ld=$2 ; shift 2 ;;
    --cov_names) cov_names=$2 ; shift 2 ;;
    --pheno_name) pheno_name=$2 ; shift 2 ;;
    --by_chr_sumstats) by_chr_sumstats=$2 ; shift 2 ;;
    --pop1) pop1=$2 ; shift 2 ;;
    --pop2) pop2=$2 ; shift 2 ;;
    --thinned_snplist) thinned_snplist=$2 ; shift 2 ;;
    --n_max_locus) n_max_locus=$2 ; shift 2 ;;
    --N_pop1) N_pop1=$2 ; shift 2 ;;
    --N_pop2) N_pop2=$2 ; shift 2 ;;
    --prop_train) prop_train=$2 ; shift 2 ;;
    --prop_test) prop_test=$2 ; shift 2 ;;
    --ranking) ranking=$2 ; shift 2 ;;
    --ids_col) ids_col=$2 ; shift 2 ;;
    --sumstats_p) sumstats_p=$2 ; shift 2 ;;
    --sumstats_snpID) sumstats_snpID=$2 ; shift 2 ;;
    --sumstats_beta) sumstats_beta=$2 ; shift 2 ;;
    --sumstats_allele1) sumstats_allele1=$2 ; shift 2 ;;
    --sumstats_allele0) sumstats_allele0=$2 ; shift 2 ;;
    --sumstats_n) sumstats_n=$2 ; shift 2 ;;
    --sumstats_se) sumstats_se=$2 ; shift 2 ;;
    --sumstats_frq) sumstats_frq=$2 ; shift 2 ;;
    --strand_check) strand_check=$2 ; shift 2 ;;
    --n_folds) n_folds=$2 ; shift 2 ;;
    --clean) clean=$2 ; shift 2 ;;
    --fst) fst=$2 ; shift 2 ;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    # If invalid options were passed, then getopt should have reported an error,
    # which we checked as VALID_ARGUMENTS when getopt was called...
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

if [ $by_chr_ld = 0 ]
then
    by_chr_ld=$by_chr
fi
if [ $by_chr_target = 0 ]
then
    by_chr_target=$by_chr
fi
if [ $blockdir = 0 ]
then
   blockdir=$outdir/blocks
fi

echo "Options in effect:"
echo "outdir  : $outdir"
echo "n_cores : $n_cores"
echo "pop1_ld_ids : $pop1_ld_ids"
echo "pop2_ld_ids : $pop2_ld_ids"
echo "pop1_ld_bfile : $pop1_ld_bfile"
echo "pop2_ld_bfile : $pop2_ld_bfile"
echo "pop1_sumstats : $pop1_sumstats"
echo "pop2_sumstats : $pop2_sumstats"
echo "outdir : $outdir"
echo "pop1_qc_snplist : $pop1_qc_snplist"
echo "pop2_qc_snplist : $pop2_qc_snplist"
echo "do_block_pop1 : $do_block_pop1"
echo "do_block_pop2 : $do_block_pop2"
echo "do_sumstat_pop1 : $do_sumstat_pop1"
echo "do_sumstat_pop2 : $do_sumstat_pop2"
echo "do_clump_pop1 : $do_clump_pop1"
echo "do_est_beta_pop1 : $do_est_beta_pop1"
echo "do_sumstat_ensembl_pop1 : $do_sumstat_ensembl_pop1"
echo "do_est_beta_pop1_precision : $do_est_beta_pop1_precision"
echo "do_est_beta_InformPrior : $do_est_beta_InformPrior"
echo "do_clump_pop2 : $do_clump_pop2"
echo "do_est_beta_pop2 : $do_est_beta_pop2"
echo "do_sumstat_ensembl_pop2 : $do_sumstat_ensembl_pop2"
echo "by_chr_target : $by_chr_target"
echo "by_chr_ld : $by_chr_ld"
echo "by_chr_sumstats : $by_chr_sumstats"
echo "pheno_name : $pheno_name"
echo "cov_names : $cov_names"
echo "pop1 : $pop1"
echo "pop2 : $pop2"
echo "n_max_locus : $n_max_locus"
echo "thinned_snplist : $thinned_snplist"
echo "ranking : $ranking"
echo "ids_col : $ids_col"
echo "strand_check : $strand_check"
echo "fst : $fst"
echo "n_folds : $n_folds"
echo ""

mkdir -p $outdir
mkdir -p $outdir/$pop1
mkdir -p $outdir/$pop2
mkdir -p $blockdir
mkdir -p $blockdir/$pop1
mkdir -p $blockdir/$pop2

if [ $ranking != "pv" ] && [ $ranking != "pv.minP" ] && [ $ranking != "pv.ftest" ] && [ $ranking != "thinned.pv.ftest" ] && [ $ranking != "f.stat" ] && [ $ranking != "thinned.f.stat" ]
then
    echo "Invalid argument, ranking="$ranking
    exit
fi

n_iter=5

for ((iter = 1; iter <= ${n_folds}; iter++))
do
    mkdir -p $outdir/$pop1/fold$iter
    mkdir -p $outdir/$pop1/fold$iter/clump
    mkdir -p $outdir/$pop1/fold$iter/models
    mkdir -p $outdir/$pop1/fold$iter/models/lambda
    mkdir -p $outdir/$pop1/fold$iter/sumstat_subset
    
    mkdir -p $outdir/$pop2/fold$iter
    mkdir -p $outdir/$pop2/fold$iter/clump
    mkdir -p $outdir/$pop2/fold$iter/models
    mkdir -p $outdir/$pop2/fold$iter/models/lambda
    mkdir -p $outdir/$pop2/fold$iter/sumstat_subset
done

if [ $do_block_pop1 -eq 1 ]
then
    for chr in {1..22}
    do
	bfile1=$pop1_ld_bfile
	if [ $by_chr -eq 1 ]
	then
	    bfile1=$pop1_ld_bfile$chr
	fi
	safe_rm $blockdir/$pop1/chr$chr.blocks*
	plink --bfile $bfile1 \
	      --chr $chr \
	      --keep $pop1_ld_ids \
	      --extract $pop1_qc_snplist \
	      --maf 0.001 \
	      --blocks 'no-pheno-req' 'no-small-max-span' \
	      --blocks-min-maf 0.001 \
	      --blocks-strong-lowci 0.51 \
	      --blocks-strong-highci 0.831 \
	      --blocks-recomb-highci 0.55 \
	      --blocks-max-kb 500 \
	      --blocks-inform-frac 0.1 \
	      --out $blockdir/${pop1}/chr${chr}
	gzip $blockdir/$pop1/chr$chr.blocks.det
    done
fi

if [ $do_block_pop2 -eq 1 ]
then
    for chr in {1..22}
    do
	bfile1=$pop2_ld_bfile
	if [ $by_chr -eq 1 ]
	then
	    bfile1=$pop2_ld_bfile$chr
	fi
	safe_rm $blockdir/${pop2}/chr$chr.blocks*
	plink --bfile $bfile1 \
	      --chr $chr \
	      --keep $pop2_ld_ids \
	      --extract $pop2_qc_snplist \
	      --maf 0.001 \
	      --blocks 'no-pheno-req' 'no-small-max-span' \
	      --blocks-min-maf 0.001 \
	      --blocks-strong-lowci 0.51 \
	      --blocks-strong-highci 0.831 \
	      --blocks-recomb-highci 0.55 \
	      --blocks-max-kb 1000 \
	      --blocks-inform-frac 0.1 \
	      --out $blockdir/${pop2}/chr${chr}
	gzip $blockdir/${pop2}/chr$chr.blocks.det
    done
fi

if [ $do_sumstat_pop1 -eq 1 ]
then
    for ((iter = 1; iter <= ${n_folds}; iter++))
    do
	safe_rm -Rf $outdir/$pop1/fold$iter/sumstat_subset/*
    done
    Rscript --vanilla $RSCRIPTS"/"make_sumstats_subset.R \
 	    --fpath $FPATH \
	    --blockdir $blockdir/$pop1 \
	    --workdir $outdir/$pop1 \
	    --sumstats $pop1_sumstats \
	    --ld.ids $pop1_ld_ids \
	    --bfile $pop1_ld_bfile \
	    --sumstats.snpID $sumstats_snpID \
	    --sumstats.betaID $sumstats_beta \
 	    --sumstats.allele1ID $sumstats_allele1 \
	    --sumstats.allele0ID $sumstats_allele0 \
	    --sumstats.P $sumstats_p \
	    --qc.snplist $pop1_qc_snplist \
	    --N.pop $N_pop1 \
	    --prop.train $prop_train \
	    --prop.test $prop_test \
	    --n.cores $n_cores \
	    --by.chr.sumstats $by_chr_sumstats \
	    --strand.check $strand_check \
	    --n.folds $n_folds \
	    --by.chr $by_chr_ld
fi

if [ $do_sumstat_pop2 -eq 1 ]
then
    for ((iter = 1; iter <= ${n_folds}; iter++))
    do
	safe_rm -Rf $outdir/$pop2/fold$iter/sumstat_subset/*
    done
    Rscript --vanilla $RSCRIPTS"/"make_sumstats_subset.R \
 	    --fpath $FPATH \
	    --blockdir $blockdir/$pop2 \
	    --workdir $outdir/$pop2 \
	    --sumstats $pop2_sumstats \
	    --ld.ids $pop2_ld_ids \
	    --bfile $pop2_ld_bfile \
	    --sumstats.snpID $sumstats_snpID \
	    --sumstats.betaID $sumstats_beta \
 	    --sumstats.allele1ID $sumstats_allele1 \
	    --sumstats.allele0ID $sumstats_allele0 \
	    --sumstats.P $sumstats_p \
	    --qc.snplist $pop2_qc_snplist \
	    --N.pop $N_pop2 \
	    --prop.train $prop_train \
	    --prop.test $prop_test \
	    --n.cores $n_cores \
	    --by.chr.sumstats $by_chr_sumstats \
	    --strand.check $strand_check \
	    --n.folds $n_folds \
	    --by.chr $by_chr_ld
fi

for ((iter = 1; iter <= ${n_folds}; iter++))
do

if [ $do_clump_pop1 -eq 1 ]
then
    for chr in {1..22}
    do
	bfile1=$pop1_ld_bfile$chr
	pop1_sumstats1=$outdir/$pop1/fold$iter/sumstat_subset/chr${chr}.dat.gz
	plink --bfile $bfile1 \
	      --chr $chr \
	      --clump $pop1_sumstats1 \
	      --clump-field P --clump-snp-field SNP \
	      --clump-p1 1e-1 --clump-p2 1e-1 --clump-kb 1000 --clump-r2 0.01 \
	      --keep $pop1_ld_ids \
	      --extract $pop1_qc_snplist \
	      --maf 0.001 \
	      --out $outdir/${pop1}/fold$iter/clump/_${chr}
	safe_rm -f $outdir/${pop1}/fold$iter/clump/_${chr}.clumped.gz
	gzip $outdir/${pop1}/fold$iter/clump/_${chr}.clumped
    done
fi

if [ $do_est_beta_pop1 -eq 1  ]
then
    safe_rm -f $outdir/$pop1/fold$iter/models/stage1*
    Rscript --vanilla $RSCRIPTS"/"est_beta_bychr.R \
 	    --fpath $FPATH \
	    --clump.stem $outdir/$pop1/fold$iter/clump/ \
	    --sumstats $outdir/$pop1/fold$iter/sumstat_subset/chr \
	    --by.chr.sumstats .dat.gz \
	    --thinned.snplist $thinned_snplist \
	    --n.max.locus $n_max_locus \
	    --ld.ids $pop1_ld_ids \
	    --beta.stem $outdir/$pop1/fold$iter/models/stage1 \
	    --bfile $pop1_ld_bfile \
	    --lambda 0.1,0.2,0.5,1,2,5,10,20 \
	    --S 0,0.25,0.5,0.75,1 \
	    --sumstats.snpID SNP \
	    --sumstats.betaID BETA \
	    --sumstats.allele1ID ALLELE1 \
	    --sumstats.allele0ID ALLELE0 \
	    --n.cores $n_cores \
	    --by.chr $by_chr_ld \
	    --strand.check $strand_check
fi

if [ $do_sumstat_ensembl_pop1 -eq 1  ]
then
    safe_rm -f $outdir/${pop1}/fold${iter}/*.dat
    Rscript --vanilla $RSCRIPTS"/"all_snp_weights.R \
 	    --fpath $FPATH \
	    --stage1  $outdir/${pop1}/fold${iter}/models/stage1 \
	    --fold $iter \
	    --blockdir $blockdir/$pop1 \
	    --workdir $outdir/${pop1} \
	    --bfile $pop1_ld_bfile \
	    --ld.ids $pop1_ld_ids \
	    --strand.check $strand_check \
	    --N.pop $N_pop1 \
	    --prop.train $prop_train \
	    --prop.test $prop_test \
	    --n.cores $n_cores

    if [ $clean -eq 1 ]
    then
	safe_rm $outdir/$pop1/fold${iter}/models/stage1*
    fi
fi

if [ $do_est_beta_pop1_precision -eq 1  ]
then
    safe_rm -f $outdir/${pop1}/fold${iter}/models/prior*
    safe_rm -f $outdir/${pop1}/fold${iter}/models/lambda/rs*.gz
    Rscript --vanilla $RSCRIPTS"/"est_beta_bychr.R \
 	    --fpath $FPATH \
	    --clump.stem $outdir/$pop1/fold$iter/clump/ \
	    --sumstats $outdir/$pop1/fold$iter/sumstat_subset/chr \
	    --by.chr.sumstats .dat.gz \
	    --thinned.snplist $thinned_snplist \
	    --n.max.locus $n_max_locus \
	    --ld.ids $pop1_ld_ids \
	    --bfile $pop1_ld_bfile \
	    --sumstats.snpID SNP \
	    --sumstats.betaID BETA \
	    --sumstats.allele1ID ALLELE1 \
	    --sumstats.allele0ID ALLELE0 \
	    --n.cores $n_cores \
	    --by.chr $by_chr_ld \
	    --strand.check $strand_check \
	    --beta.stem $outdir/${pop1}/fold${iter}/models/prior \
	    --param.file $outdir/${pop1}/fold$iter/best_model_params.dat \
	    --precision TRUE

    if [ $clean -eq 1 ]
    then
	safe_rm $outdir/$pop1/fold$iter/sumstat_subset/chr*
	safe_rm $outdir/$pop1/fold$iter/clump/_*
    fi
fi

if [ $do_est_beta_InformPrior -eq 1  ]
then
    safe_rm -f $outdir/${pop2}/fold$iter/models/stage2*
    Rscript --vanilla $RSCRIPTS"/"est_beta_InformPrior_bychr.R \
 	    --fpath $FPATH \
	    --sumstats $outdir/$pop2/fold$iter/sumstat_subset/chr  \
	    --ld.ids $pop2_ld_ids \
	    --prior $outdir/${pop1}/fold${iter}/models/prior \
	    --param.file $outdir/${pop1}/fold${iter}/best_model_params.dat \
	    --beta.stem $outdir/${pop2}/fold$iter/models/stage2 \
	    --bfile $pop2_ld_bfile \
	    --fst $fst \
	    --sumstats.snpID SNP \
	    --sumstats.betaID BETA \
	    --sumstats.allele1ID ALLELE1 \
	    --sumstats.allele0ID ALLELE0 \
	    --sumstats.P P \
	    --N.pop1 $N_pop1 \
	    --N.pop2 $N_pop2 \
	    --n.cores $n_cores \
	    --by.chr.sumstats .dat.gz \
	    --ranking pv \
	    --strand.check $strand_check \
	    --by.chr $by_chr_ld

    if [ $clean -eq 1 ]
    then
	safe_rm -Rf $outdir/$pop1/fold${iter}/models/lambda
	safe_rm $outdir/$pop1/fold${iter}/models/prior*
    fi
fi

if [ $do_clump_pop2 -eq 1 ]
then
    for chr in {1..22}
    do
	bfile1=$pop2_ld_bfile$chr
	pop2_sumstats1=$outdir/$pop2/fold$iter/sumstat_subset/chr${chr}.dat.gz
	plink --bfile $bfile1 \
	      --chr $chr \
	      --clump $pop2_sumstats1 \
	      --clump-field P --clump-snp-field SNP \
	      --clump-p1 1e-1 --clump-p2 1e-1 --clump-kb 1000 --clump-r2 0.01 \
	      --keep $pop2_ld_ids \
	      --extract $pop2_qc_snplist \
	      --maf 0.001 \
	      --out $outdir/${pop2}/fold$iter/clump/_${chr}
	safe_rm -f $outdir/${pop2}/fold$iter/clump/_${chr}.clumped.gz
	gzip $outdir/${pop2}/fold$iter/clump/_${chr}.clumped
    done
fi

if [ $do_est_beta_pop2 -eq 1  ]
then
    safe_rm -f $outdir/$pop2/fold$iter/models/stage1*
    Rscript --vanilla $RSCRIPTS"/"est_beta_bychr.R \
 	    --fpath $FPATH \
	    --clump.stem $outdir/$pop2/fold$iter/clump/ \
	    --sumstats $outdir/$pop2/fold$iter/sumstat_subset/chr \
	    --thinned.snplist $thinned_snplist \
	    --n.max.locus $n_max_locus \
	    --ld.ids $pop2_ld_ids \
	    --beta.stem $outdir/$pop2/fold$iter/models/stage1 \
	    --bfile $pop2_ld_bfile \
	    --lambda 0.2,0.5,1,2,5,10,20,50 \
	    --S 0,0.25,0.5,0.75,1 \
	    --sumstats.snpID SNP \
	    --sumstats.betaID BETA \
	    --sumstats.allele1ID ALLELE1 \
	    --sumstats.allele0ID ALLELE0 \
	    --n.cores $n_cores \
	    --by.chr $by_chr_ld \
	    --by.chr.sumstats .dat.gz \
	    --strand.check $strand_check

    if [ $clean -eq 1 ]
    then
	safe_rm $outdir/$pop2/fold$iter/clump/_*
    fi
fi

if [ $do_sumstat_ensembl_pop2 -eq 1  ]
then
    safe_rm -f $outdir/$pop2/fold$iter/*.dat
    Rscript --vanilla $RSCRIPTS"/"all_snp_weights.R \
 	    --fpath $FPATH \
	    --stage1  $outdir/$pop2/fold$iter/models/stage1 \
	    --stage2 $outdir/$pop2/fold$iter/models/stage2 \
	    --blockdir $blockdir/$pop2 \
	    --workdir $outdir/$pop2 \
	    --bfile $pop2_ld_bfile \
	    --ld.ids $pop2_ld_ids \
	    --strand.check $strand_check \
	    --N.pop $N_pop2 \
	    --prop.train $prop_train \
	    --prop.test $prop_test \
	    --fold $iter \
	    --n.cores $n_cores

    if [ $clean -eq 1 ]
    then
	safe_rm $outdir/$pop2/fold$iter/sumstat_subset/chr*
	safe_rm $outdir/$pop2/fold${iter}/models/stage*
    fi
fi
done

if [ $do_pool -eq 1  ]
then
    safe_rm $outdir/${pop2}/snp_weights_*_model.dat
    Rscript --vanilla $RSCRIPTS"/"pool_snp_scores.R \
	    --n.folds $n_folds \
	    --workdir $outdir/$pop2
fi

#!/bin/bash
# update key
# Set some default values:

RSCRIPTS=~/BridgePRS/src/Rscripts
FPATH=$RSCRIPTS"/functions.R"

pheno_name="y"
ext=""
n_cores=1
do_clump_pop1=1
do_est_beta_pop1=1
do_predict_pop1=1
do_est_beta_pop1_precision=1
do_est_beta_InformPrior=1
do_predict_pop2_stage2=1
do_clump_pop2=1
do_est_beta_pop2=1
do_predict_pop2=1
do_combine=1
by_chr=0
by_chr_ld=0
by_chr_target=0
by_chr_sumstats=0
indir=0
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
binary=0

prop_train=0.7
prop_test=0.15

#RSCRIPTS=~/BridgePRS/src/Rscripts
#FPATH=$RSCRIPTS"/functions.R"


#usage()
#{
#  echo "Usage: ridgePRS [ -n | --n_cores n_cores ] [ -b | --bfile bfile ]
#                        [ -o | --outdir outdir ] 
#                        [ -v | --eur_ld_ids eur_ld_ids ] filename(s)"
#  exit 2
#}

PARSED_ARGUMENTS=$(getopt -a -n ridgePRS -o b:o:n:c:d:e:f:g:h:i:j:k:l:m:p:q:r:s:t:u:v:w:x:y:z:1:2:3:4:5:6:7: --long outdir:,n_cores:,pop1_ld_ids:,pop2_ld_ids:,pop1_ld_bfile:,pop2_ld_bfile:,pop1_sumstats:,pop2_sumstats:,pop1_valid_data:,pop2_valid_data:,pop1_test_data:,pop2_test_data:,pop1_qc_snplist:,pop2_qc_snplist:,do_clump_pop1:,do_est_beta_pop1:,,do_est_beta_pop1_precision:,do_est_beta_InformPrior:,do_clump_pop2:,do_est_beta_pop2:,do_est_beta_pop2:,do_block_pop1:,do_sumstat_pop1:,do_block_pop2:,do_sumstat_pop2:,do_sumstat_ensembl_pop1:,do_sumstat_ensembl_pop2:,by_chr:,by_chr_target:,by_chr_ld:,cov_names:,pheno_name:,indir:,by_chr_sumstats:,pop2:,thinned_snplist:,n_max_locus:,ranking:,pop1:,ids_col:,sumstats_snpID:,sumstats_beta:,sumstats_allele1:,sumstats_allele0:,sumstats_p:,sumstats_n:,sumstats_se:,sumstats_frq:,strand_check:,fst:,binary:,N_pop1:,N_pop2: -- "$@")
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
    --by_chr) by_chr=$2 ; shift 2 ;;
    --by_chr_target) by_chr_target=$2 ; shift 2 ;;
    --by_chr_ld) by_chr_ld=$2 ; shift 2 ;;
    --cov_names) cov_names=$2 ; shift 2 ;;
    --pheno_name) pheno_name=$2 ; shift 2 ;;
    --indir) indir=$2 ; shift 2 ;;
    --by_chr_sumstats) by_chr_sumstats=$2 ; shift 2 ;;
    --pop1) pop1=$2 ; shift 2 ;;
    --pop2) pop2=$2 ; shift 2 ;;
    --thinned_snplist) thinned_snplist=$2 ; shift 2 ;;
    --n_max_locus) n_max_locus=$2 ; shift 2 ;;
    --N_pop1) N_pop1=$2 ; shift 2 ;;
    --N_pop2) N_pop2=$2 ; shift 2 ;;
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
    --fst) fst=$2 ; shift 2 ;;
    --binary) binary=$2 ; shift 2 ;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    # If invalid options were passed, then getopt should have reported an error,
    # which we checked as VALID_ARGUMENTS when getopt was called...
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

if [ $indir = 0 ]
then
    indir=$outdir
fi
if [ $by_chr_ld = 0 ]
then
    by_chr_ld=$by_chr
fi
if [ $by_chr_target = 0 ]
then
    by_chr_target=$by_chr
fi
if [ $binary_pop1 = 0 ]
then
    binary_pop1=$binary
fi
if [ $binary_pop2 = 0 ]
then
    binary_pop2=$binary
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
echo "pop1_valid_data : $pop1_valid_data"
echo "pop2_valid_data : $pop2_valid_data"
echo "pop1_test_data : $pop1_test_data"
echo "pop2_test_data : $pop2_test_data"
echo "outdir : $outdir"
echo "pop1_qc_snplist : $pop1_qc_snplist"
echo "pop2_qc_snplist : $pop2_qc_snplist"
echo "do_clump_pop1 : $do_clump_pop1"
echo "do_est_beta_pop1 : $do_est_beta_pop1"
echo "do_est_beta_pop1_precision : $do_est_beta_pop1_precision"
echo "do_est_beta_InformPrior : $do_est_beta_InformPrior"
echo "do_clump_pop2 : $do_clump_pop2"
echo "do_est_beta_pop2 : $do_est_beta_pop2"
echo "do_combine : $do_combine"
echo "by_chr_target : $by_chr_target"
echo "by_chr_ld : $by_chr_ld"
echo "by_chr_sumstats : $by_chr_sumstats"
echo "pheno_name : $pheno_name"
echo "cov_names : $cov_names"
echo "indir : $indir"
echo "pop1 : $pop1"
echo "pop2 : $pop2"
echo "n_max_locus : $n_max_locus"
echo "thinned_snplist : $thinned_snplist"
echo "ranking : $ranking"
echo "ids_col : $ids_col"
echo "strand_check : $strand_check"
echo "fst : $fst"
echo "binary : $binary"
echo ""

mkdir -p $outdir
mkdir -p $outdir/$pop1/blocks
mkdir -p $outdir/$pop2/blocks
mkdir -p $outdir/$pop1/sumstat_subset
mkdir -p $outdir/$pop2/sumstat_subset
mkdir -p $outdir/clump
mkdir -p $outdir/models
mkdir -p $outdir/models/lambda

if [ $ranking != "pv" ] && [ $ranking != "pv.minP" ] && [ $ranking != "pv.ftest" ] && [ $ranking != "thinned.pv.ftest" ] && [ $ranking != "f.stat" ] && [ $ranking != "thinned.f.stat" ]
then
    echo "Invalid argument, ranking="$ranking
    exit
fi

if [ $indir == "" ]
then
    indir=$outdir
fi

if [ $do_block_pop1 -eq 1 ]
then
    for chr in {1..22}
    do
	bfile1=$pop1_ld_bfile
	if [ $by_chr -eq 1 ]
	then
	    bfile1=$pop1_ld_bfile$chr
	fi
	pop1_sumstats1=$pop1_sumstats
	if [ $by_chr_sumstats != 0 ]
	then
	    pop1_sumstats1=$pop1_sumstats$chr$by_chr_sumstats
	fi
	rm -f $outdir/$pop1/blocks/chr$chr.blocks*
	col_num=$(zcat $pop1_sumstats1 | 
		      head -1 | tr ' ' '\n' | cat -n | grep ID | awk '{print $1}')
	zcat $pop1_sumstats1 | cut -d " " -f $col_num \
				   > $outdir/$pop1/blocks/chr${chr}.snplist
	plink --bfile $bfile1 \
	      --chr $chr \
	      --keep $pop1_ld_ids \
	      --extract $outdir/$pop1/blocks/chr${chr}.snplist \
	      --maf 0.001 \
	      --blocks 'no-pheno-req' 'no-small-max-span' \
	      --blocks-min-maf 0.00001 \
	      --blocks-strong-lowci 0.51 \
	      --blocks-strong-highci 0.831 \
	      --blocks-recomb-highci 0.55 \
	      --blocks-max-kb 1000 \
	      --blocks-inform-frac 0.1 \
	      --out $outdir/${pop1}/blocks/chr${chr}
	rm $outdir/$pop1/blocks/chr$chr.blocks
	gzip $outdir/$pop1/blocks/chr$chr.blocks.det
    done
fi

if [ $do_sumstat_pop1 -eq 1 ]
then
    rm $outdir/$pop1/sumstat_subset/*
    Rscript --vanilla $RSCRIPTS"/"make_sumstats_subset.R \
 	    --fpath $FPATH \
	    --workdir $outdir/$pop1 \
	    --sumstats $pop1_sumstats \
	    --ld.ids $pop1_ld_ids \
	    --bfile $pop1_ld_bfile \
	    --sumstats.snpID $sumstats_snpID \
	    --sumstats.betaID $sumstats_beta \
 	    --sumstats.allele1ID $sumstats_allele1 \
	    --sumstats.allele0ID $sumstats_allele0 \
	    --sumstats.P $sumstats_p \
	    --N.pop $N_pop1 \
	    --prop.train $prop_train \
	    --prop.test $prop_test \
	    --n.cores $n_cores \
	    --by.chr.sumstats $by_chr_sumstats \
	    --strand.check $strand_check \
	    --by.chr $by_chr_ld
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
	pop2_sumstats1=$pop2_sumstats
	if [ $by_chr_sumstats != 0 ]
	then
	    pop2_sumstats1=$pop2_sumstats$chr$by_chr_sumstats
	fi
	rm -f $outdir/$pop2/blocks/chr$chr.blocks*
	col_num=$(zcat $pop2_sumstats1 | 
		      head -1 | tr ' ' '\n' | cat -n | grep ID | awk '{print $1}')
	zcat $pop2_sumstats1 | cut -d " " -f $col_num \
				   > $outdir/${pop2}/blocks/chr${chr}.snplist
	plink --bfile $bfile1 \
	      --chr $chr \
	      --keep $pop2_ld_ids \
	      --extract $outdir/${pop2}/blocks/chr${chr}.snplist \
	      --maf 0.001 \
	      --blocks 'no-pheno-req' 'no-small-max-span' \
	      --blocks-min-maf 0.00001 \
	      --blocks-strong-lowci 0.51 \
	      --blocks-strong-highci 0.831 \
	      --blocks-recomb-highci 0.55 \
	      --blocks-max-kb 1000 \
	      --blocks-inform-frac 0.1 \
	      --out $outdir/${pop2}/blocks/chr${chr}
	rm $outdir/${pop2}/blocks/chr$chr.blocks
	gzip $outdir/${pop2}/blocks/chr$chr.blocks.det
    done
fi

if [ $do_sumstat_pop2 -eq 1 ]
then
    rm $outdir/$pop2/sumstat_subset/*
    Rscript --vanilla $RSCRIPTS"/"make_sumstats_subset.R \
 	    --fpath $FPATH \
	    --workdir $outdir/$pop2 \
	    --sumstats $pop2_sumstats \
	    --ld.ids $pop2_ld_ids \
	    --bfile $pop2_ld_bfile \
	    --sumstats.snpID $sumstats_snpID \
	    --sumstats.betaID $sumstats_beta \
 	    --sumstats.allele1ID $sumstats_allele1 \
	    --sumstats.allele0ID $sumstats_allele0 \
	    --sumstats.P $sumstats_p \
	    --N.pop $N_pop2 \
	    --prop.train $prop_train \
	    --prop.test $prop_test \
	    --n.cores $n_cores \
	    --by.chr.sumstats $by_chr_sumstats \
	    --strand.check $strand_check \
	    --by.chr $by_chr_ld
fi

if [ $do_clump_pop1 -eq 1 ]
then
    for chr in {1..22}
    do
	bfile1=$pop1_ld_bfile
	if [ $by_chr -eq 1 ]
	then
	    bfile1=$pop1_ld_bfile$chr
	fi
	pop1_sumstats1=$outdir/$pop1/sumstat_subset/chr
	if [ $by_chr_sumstats != 0 ]
	then
	    pop1_sumstats1=$outdir/$pop1/sumstat_subset/chr$chr.dat.gz
	fi
	plink --bfile $bfile1 \
	      --chr $chr \
	      --clump $pop1_sumstats1 \
	      --clump-field P --clump-snp-field SNP \
	      --clump-p1 1e-1 --clump-p2 1e-1 --clump-kb 1000 --clump-r2 0.01 \
	      --keep $pop1_ld_ids \
	      --extract $pop1_qc_snplist \
	      --maf 0.001 \
	      --out $outdir/clump/$pop1\_$chr
	rm -f $outdir/clump/$pop1\_$chr.clumped.gz
	gzip $outdir/clump/$pop1\_$chr.clumped
    done
fi

if [ $do_est_beta_pop1 -eq 1  ]
then
    rm -f $outdir/models/$pop1*
    Rscript --vanilla $RSCRIPTS"/"est_beta_bychr.R \
 	    --fpath $FPATH \
	    --clump.stem $indir/clump/$pop1 \
	    --sumstats $outdir/$pop1/sumstat_subset/chr \
	    --thinned.snplist $thinned_snplist \
	    --n.max.locus $n_max_locus \
	    --ld.ids $pop1_ld_ids \
	    --beta.stem $outdir/models/$pop1\_stage1 \
	    --bfile $pop1_ld_bfile \
	    --lambda 0.1,0.2,0.5,1,2,5,10,20 \
	    --S 0,0.25,0.5,0.75,1 \
	    --sumstats.snpID SNP \
	    --sumstats.betaID BETA \
	    --sumstats.allele1ID ALLELE1 \
	    --sumstats.allele0ID ALLELE0 \
	    --n.cores $n_cores \
	    --by.chr $by_chr_ld \
	    --by.chr.sumstats .dat.gz \
	    --strand.check $strand_check
fi

if [ $do_sumstat_ensembl_pop1 -eq 1  ]
then
    rm -f $outdir/models/$pop2\_stage1*
    Rscript --vanilla $RSCRIPTS"/"all_snp_weights.R \
 	    --fpath $FPATH \
	    --stage1  $outdir/models/$pop1\_stage1 \
	    --workdir $outdir/$pop1 \
	    --bfile $pop1_ld_bfile \
	    --ld.ids $pop1_ld_ids \
	    --strand.check $strand_check \
	    --N.pop $N_pop1 \
	    --prop.train $prop_train \
	    --prop.test $prop_test \
	    --n.cores $n_cores
fi

if [ $do_est_beta_pop1_precision -eq 1  ]
then
    rm -f $outdir/models/stage1*
    rm -f $outdir/models/lambda/rs*.gz
    Rscript --vanilla $RSCRIPTS"/"est_beta_bychr.R \
 	    --fpath $FPATH \
	    --clump.stem $indir/clump/$pop1 \
	    --sumstats $outdir/$pop1/sumstat_subset/chr \
	    --thinned.snplist $thinned_snplist \
	    --n.max.locus $n_max_locus \
	    --ld.ids $pop1_ld_ids \
	    --bfile $pop1_ld_bfile \
	    --beta.stem $outdir/models/stage1 \
	    --param.file $outdir/$pop1/best_model_params.dat \
	    --precision TRUE \
	    --sumstats.snpID SNP \
	    --sumstats.betaID BETA \
	    --sumstats.allele1ID ALLELE1 \
	    --sumstats.allele0ID ALLELE0 \
	    --n.cores $n_cores \
	    --by.chr.sumstats .dat.gz \
	    --strand.check $strand_check \
	    --by.chr $by_chr_ld
fi

if [ $do_est_beta_InformPrior -eq 1  ]
then
    rm -f $outdir/models/${pop2}_stage2*
    Rscript --vanilla $RSCRIPTS"/"est_beta_InformPrior_bychr.R \
 	    --fpath $FPATH \
	    --sumstats $outdir/$pop2/sumstat_subset/chr  \
	    --ld.ids $pop2_ld_ids \
	    --prior $outdir/models/stage1 \
	    --param.file $outdir/${pop1}/best_model_params.dat \
	    --beta.stem $outdir/models/${pop2}_stage2 \
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
fi

if [ $do_clump_pop2 -eq 1  ]
then
    for chr in {1..22}
    do
	bfile1=$pop2_ld_bfile
	if [ $by_chr -eq 1 ]
	then
	    bfile1=$pop2_ld_bfile$chr
	fi
	pop2_sumstats1=$outdir/$pop2/sumstat_subset/chr
	if [ $by_chr_sumstats != 0 ]
	then
	    pop2_sumstats1=$outdir/$pop2/sumstat_subset/chr$chr.dat.gz
	fi
	plink --bfile $bfile1 \
	      --chr $chr \
	      --clump $pop2_sumstats1 \
	      --clump-field P --clump-snp-field SNP \
	      --clump-p1 1e-2 --clump-p2 1e-1 --clump-kb 1000 --clump-r2 0.01 \
	      --keep $pop2_ld_ids \
	      --extract $pop2_qc_snplist \
	      --maf 0.001 \
	      --out $outdir/clump/$pop2\_$chr
	rm -f $outdir/clump/$pop2\_$chr.clumped.gz
	gzip $outdir/clump/$pop2\_$chr.clumped
    done
fi

if [ $do_est_beta_pop2 -eq 1  ]
then
    rm -f $outdir/models/$pop2\_stage1*
    Rscript --vanilla $RSCRIPTS"/"est_beta_bychr.R \
 	    --fpath $FPATH \
	    --clump.stem $indir/clump/$pop2 \
	    --sumstats $outdir/$pop2/sumstat_subset/chr \
	    --thinned.snplist $thinned_snplist \
	    --n.max.locus $n_max_locus \
	    --ld.ids $pop2_ld_ids \
	    --beta.stem $outdir/models/$pop2\_stage1 \
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
fi

if [ $do_sumstat_ensembl_pop2 -eq 1  ]
then
    Rscript --vanilla $RSCRIPTS"/"all_snp_weights.R \
 	    --fpath $FPATH \
	    --stage1  $outdir/models/$pop2\_stage1 \
	    --stage2 $outdir/models/$pop2\_stage2 \
	    --workdir $outdir/$pop2 \
	    --bfile $pop2_ld_bfile \
	    --ld.ids $pop2_ld_ids \
	    --strand.check $strand_check \
	    --N.pop $N_pop1 \
	    --prop.train $prop_train \
	    --prop.test $prop_test \
	    --n.cores $n_cores
fi

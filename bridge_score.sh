dir1=~/BridgePRS/out
dir2=~/BridgePRS/data
model=(best weighted)
pop=EAS

for i in {0..1}
do
    for chr in {1..22}
    do
	/hpc/users/hoggac01/plink2 --bfile $dir2/pop_${pop}/genotypes/chr$chr \
				   --score $dir1/${pop}_snp_weights_${model[$i]}_model_${iter}.dat 1 2 4 \
				   --keep $dir2/pop_${pop}/phenotypes/${pop}_valid.dat \
				   --out $dir1/preds/${pop}_${model[$i]}_${iter}_chr$chr
    done
done

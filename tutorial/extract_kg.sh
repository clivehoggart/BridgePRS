for chr in {1..22}
do
    ~/plink2 --bfile /sc/arion/projects/data-ark/1000G/phase3/PLINK/chr$chr \
	     --extract /sc/arion/projects/psychgen/projects/BridgePRS/new_data/csxsnp \
	     --make-bed \
	     --out /sc/arion/projects/psychgen/projects/BridgePRS/new_data/bfiles/chr$chr
done
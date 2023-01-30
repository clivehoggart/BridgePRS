pheno=(50 21001 30610 30220 30710 30150 30780 30130 30140 30080 30870 30860)

for i in {0..11}
do
    cov_names=Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21,PC22,PC23,PC24,PC25,PC26,PC27,PC28,PC29,PC30,PC31,PC32,PC33,PC34,PC35,PC36,PC37,PC38,PC39,PC40
    if [ ${pheno[$i]} = 30780 ] || [ ${pheno[$i]} = 30610 ] || [ ${pheno[$i]} = 30710 ] || [ ${pheno[$i]} = 30870 ] || [ ${pheno[$i]} = 30860 ]
    then
	cov_names=$cov_names,Fasting,Dilution
    fi
    echo ${pheno[$i]}
    echo $cov_names
done

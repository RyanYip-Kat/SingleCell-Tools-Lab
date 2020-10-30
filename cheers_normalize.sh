src="/home/ye/Work/BioAligment/SNP/CHEERS/CHEERS_normalize.py"
python="/home/ye/anaconda3/envs/ldsc/bin/python"

LD_DIR="/home/ye/Data/SNP_List/cheers"
FeatureCount_DIR="/home/ye/Work/BioAligment/SNP/featureCount"

OUT="./LD_Cheers_Normalize/"
if [ ! -d $OUT ];
then
        mkdir $OUT
fi

$python $src --input $FeatureCount_DIR  --outdir $OUT --prefix "GDT_Cell" 

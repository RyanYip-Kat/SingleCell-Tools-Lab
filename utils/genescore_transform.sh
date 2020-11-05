transform_src="/home/ye/Work/R/scATAC/ArchR/utils/transform_h5.sh"
gscore="/home/ye/Work/R/scATAC/ArchR/utils/GeneScoreMatrix2H5.R"
rscript="/home/ye/anaconda3/envs/scatac/bin/Rscript"

archr=$1
out=$2
CRH5=$3
if [ ! -d $out ]
then
	mkdir -p $out
fi

echo "### ArchR Project From : $archr,and Write 10x Matrix into : $out"
$rscript $gscore --project $archr --outdir $out
cd $out/GeneScoreMatrix  && awk '{print $1"\t"$2"\tGene Expression"}' genes.tsv > features.tsv && rm genes.tsv && sed -i 's/\./\-/' barcodes.tsv && gzip *
cd -

echo "### Transform H5 File into cellrange style"
bash $transform_src --target-matrix ${out}/GeneScoreMatrix --cellranger-matrix $CRH5 --output-dir $out

# cellranger reanalyze --id test  --matrix Cytof-H5/VKH/filtered_feature_bc_matrix_gold.h5 

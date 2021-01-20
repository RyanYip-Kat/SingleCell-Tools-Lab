transform="/home/ye/Work/R/Cytof/cellranger/H5_To_CRH5.R"
rscript="/home/ye/anaconda3/envs/r-base/bin/Rscript"
src="/home/ye/Work/R/Cytof/cellranger/process_fcs_not_same_beatname.R"
cellranger="/home/ye/Software/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger"

out=$1
id=$2
CRH5="/home/ye/Data/dataset/GT-GD/5RNA/GT/G-T-01-10X5/outs/filtered_feature_bc_matrix.h5"
if [ ! -d $out ]
then
	mkdir -p $out
fi


$rscript $src --samples sample_un.csv --config config_uns.csv --number 5000 --batch_correct --outdir $out  #  $out/matrix,$out/filtered_feature_bc_matrix.h5  
cd ${out}/matrix  && awk '{print $1"\t"$2"\tGene Expression"}' genes.tsv > features.tsv && rm genes.tsv && sed -i 's/\./\-/' barcodes.tsv && gzip *
cd -
bash $transform --target-matrix ${out}/matrix --cellranger-matrix $CRH5 --output-dir $out
$cellranger reanalyze --id $id --matrix ${out}/filtered_feature_bc_matrix_gold.h5

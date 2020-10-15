path="/home/ye/Work/Python/SingleCell/Project/CEpi_Mouse/cytof-R"
seurat=$1
outdir=$2
port=$3

echo "# Step 1 : Prepare data"
/home/ye/anaconda3/envs/r-base/bin/Rscript $path/Seurat_prepare_anndata.R --seurat $seurat --outdir $outdir

echo "# Step 2 : Convert seurat into anndata"
/home/ye/anaconda3/envs/scanpy/bin/python  $path/Seurat_to_Anndata.py --seurat_dir $outdir --outdir $outdir

#echo "# Step 3 : "
#/home/ye/anaconda3/envs/pegasus/bin/cirro launch $outdir/seurat_adata.h5ad --host 10.100.44.197 --port $port --no-open


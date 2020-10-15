data=$1
column=$2 # louvain_labels:3,6 --subset-selection Condition:CB_nonmix
subset=$3
out=$4

if [ ! -d $out ]
then
        mkdir -p $out
fi

echo "--subset-selection $column:$subset" 
/home/ye/anaconda3/envs/pegasus/bin/pegasus subcluster -p 8 \
	--correct-batch-effect \
	--batch-group-by "Status" \
	--correction-method "harmony" \
	--knn-K 30 \
        --louvain \
        --louvain-class-label "louvain" \
        --leiden \
        --leiden-class-label "leiden" \
	--tsne \
	--umap \
	--diffmap \
	--diffmap-to-3d \
	--output-loom \
        --select-hvf-ngenes 2000 \
	--output-seurat-compatible \
	--pca-robust --subset-selection "$column:$subset" $data $out/adata


#pegasus subcluster -p 4 --correct-batch-effect  --subset-selection louvain_labels:3,6  --tsne --louvain --umap --leiden --pca-robust --select-hvf-ngenes  2000 Cornea/human/cluster.h5ad  Cornea/human/subset/adata


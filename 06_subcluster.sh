data=$1
column=$2 # louvain_labels:3,6 --subset-selection Condition:CB_nonmix
subset=$3
out=$4

echo "--subset-selection $column:$subset" 
/home/ye/anaconda3/envs/pegasus/bin/pegasus subcluster -p 4 \
	--subset-selection "$column:$subset" \
	--correct-batch-effect \
	--louvain \
	--leiden \
	--tsne \
	--umap \
	--diffmap \
	--output-loom \
        --select-hvf-ngenes 2000 \
	--pca-robust $data $out


#pegasus subcluster -p 4 --correct-batch-effect  --subset-selection louvain_labels:3,6  --tsne --louvain --umap --leiden --pca-robust --select-hvf-ngenes  2000 Cornea/human/cluster.h5ad  Cornea/human/subset/adata


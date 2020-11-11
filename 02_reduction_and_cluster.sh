data=$1
batch=$2
out=$3

/home/ye/anaconda3/envs/pegasus/bin/pegasus cluster $data $out -p 4 \
	--correct-batch-effect \
	--batch-group-by $batch \
	--correction-method "harmony" \
	--louvain \
	--leiden \
	--tsne \
	--umap \
	--diffmap \
	--mito-prefix "MT-" \
	--percent-mito 6.0 \
        --select-hvf-ngenes 2000 \
	--min-genes 3 \
	--max-genes 7500 \
	--plot-hvf \
	--pca-robust \
	--output-seurat-compatible \
	--output-loom



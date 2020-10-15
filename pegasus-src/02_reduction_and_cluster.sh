data=$1
batch=$2
out=$3

# --correction-method "L/S" or "harmony"
/home/ye/anaconda3/envs/pegasus/bin/pegasus cluster $data $out -p 16 \
	--correct-batch-effect \
	--batch-group-by $batch \
	--correction-method "harmony" \
	--knn-K 30 \
	--kBET \
	--kBET-batch $batch \
	--louvain \
	--louvain-class-label "louvain" \
	--leiden \
	--leiden-class-label "leiden" \
	--tsne \
	--umap \
	--diffmap \
	--diffmap-to-3d \
	--mito-prefix "mt-" \
	--percent-mito 6 \
        --select-hvf-ngenes 2000 \
	--min-genes 2000 \
	--max-genes 8000 \
	--plot-hvf \
	--plot-filtration-results \
	--pca-robust \
	--output-seurat-compatible \
	--output-loom



mat=$1 # linear_gene_expression_matrix.tsv
out=$2 

if [ ! -d $out ]
then
        mkdir -p "$out/compass"
	mkdir -p "$out/tmp"
fi

/home/ye/anaconda3/envs/BulkBio/bin/compass --data $mat \
       	--num-processes 16 \
        --species "homo_sapiens" \
	--output-dir "$out/compass" \
	--temp-dir  "$out/tmp"  \




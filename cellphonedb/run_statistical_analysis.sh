status=$1
meta="output/cellphonedb/$status/cell_meta.txt"
counts="output/cellphonedb/$status/cell_counts.txt"

/home/ye/anaconda3/envs/CellPhoneDB/bin/cellphonedb method statistical_analysis $meta $counts \
	--project-name "out" \
	--output-path output/cellphonedb/$status \
	--threads 12 \
	--pvalue 0.05 

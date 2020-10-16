status=$1
meta="output/cellphonedb/$status/cell_meta.txt"
counts="output/cellphonedb/$status/cell_counts.txt"
db="/home/ye/Work/Python/SingleCell/Project/cellphonedb/db/cellphonedb_user_2020-07-04-13_03.db"

/home/ye/anaconda3/envs/CellPhoneDB/bin/cellphonedb method statistical_analysis $meta $counts \
	--database $db \
	--project-name "out" \
	--output-path output/cellphonedb/$status \
	--threads 12 \
	--pvalue 0.05 

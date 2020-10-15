data=$1
out=$2

/home/ye/anaconda3/envs/pegasus/bin/pegasus annotate_cluster  $data $out/annotation.txt --marker-file human_immune \
	--de-test  "t" \
	--de-alpha  0.05 \
	--de-key "de_res" \
	--minimum-report-score 0.5 \
	--do-not-use-non-de-genes 


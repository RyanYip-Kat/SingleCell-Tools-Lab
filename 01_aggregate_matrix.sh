sample_csv=$1
out=$2
/home/ye/anaconda3/envs/pegasus/bin/pegasus  aggregate_matrix --attributes Sample,Source,Platform,Donor  $sample_csv $out

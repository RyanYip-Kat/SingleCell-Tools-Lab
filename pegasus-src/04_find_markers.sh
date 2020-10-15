data=$1
label=$2  #louvain_labels ,leiden_labels
out=$3

/home/ye/anaconda3/envs/pegasus/bin/pegasus find_markers --labels $label --remove-ribo  --min-gain 1.0 -p 4 $data $out/${label}_markers.xlsx


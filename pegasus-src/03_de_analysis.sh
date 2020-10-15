data=$1
label=$2  #louvain_labels ,leiden_labels
out=$3
/home/ye/anaconda3/envs/pegasus/bin/pegasus de_analysis $data $out/${label}_de_markers.xlsx -p 4 \
	--labels $label --auc --t --fisher --mwu --ndigits 4 --alpha 0.05


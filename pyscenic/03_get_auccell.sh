#pyscenic aucell ./loom/aging-11/BC/scenic.loom ./loom/aging-11/BC/reg.csv --output ./loom/aging-11/BC/pyscenic_output.loom --num_workers 12
loom=$1
reg=$2
out=$3
/home/ye/anaconda3/envs/scanpy/bin/pyscenic aucell  $loom $reg --output $out/pyscenic_output.loom --num_workers 12 


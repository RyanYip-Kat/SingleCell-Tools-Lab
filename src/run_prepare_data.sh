outdir="metric_qc_plot"
cat ./sample_path.txt | while read id
do
	arr=(${id})
	sample=${arr[0]}
	path=${arr[1]}
	echo "$sample ------ $path"
        /home/ye/anaconda3/envs/r-base/bin/Rscript ./preprepare_data.R --path $path --outdir ${outdir}/$sample 	
done

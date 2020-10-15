outdir="output"
cat ./sample_path.txt | while read id
do
	arr=(${id})
	sample=${arr[0]}
	path=${arr[1]}
	lower=${arr[2]}
	upper=${arr[3]}
	mt=${arr[4]}
	suffix=${arr[5]}
	echo "$sample ------ $path"
        /home/ye/anaconda3/envs/r-base/bin/Rscript ./removeCells_from_qc.R --path $path --outdir ${outdir}/$sample --lower $lower --upper $upper --mt $mt --suffix $suffix
done

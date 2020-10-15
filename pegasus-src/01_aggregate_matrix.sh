sample_csv=$1
out=$2

outdir=`dirname $out`
if [ ! -d $outdir/QC ];
then
	mkdir -p $outdir/QC
fi

/home/ye/anaconda3/envs/pegasus/bin/pegasus  aggregate_matrix --attributes Sample,Source,Platform,Donor,Status  $sample_csv $out
/home/ye/anaconda3/envs/pegasus/bin/python /home/ye/Work/Python/SingleCell/Project/CEpi_Mouse/pegasus-src/caculate_qc_from_aggregate_matrix.py --h5sc $out.h5sc --outdir $outdir/QC

data=$1
out=$3
group=$2
pegasus="/home/ye/anaconda3/envs/pegasus/bin/pegasus"

$pegasus plot composition --cluster-labels $group --attribute Donor --style normalized --not-stacked $data $out/composition.pdf
$pegasus plot scatter --basis tsne --attributes $group,Donor $data $out/scatter.pdf

$pegasus plot scatter_groups --basis tsne --cluster-labels $group --group Donor $data $out/scatter_groups.pdf
$pegasus plot scatter_genes --basis tsne --genes CD8A,CD4,CD3G,NCAM1,CD14,ITGAX,IL3RA,CD38,CD34,PPBP $data $out/genes.pdf

$pegasus plot scatter_gene_groups --basis tsne --gene CD8A --group Donor $data $out/gene_groups.pdf
$pegasus plot heatmap --cluster-labels $group --genes CD8A,CD4,CD3G,NCAM1,CD14,ITGAX,IL3RA,CD38,CD34,PPBP --show-zscore --heatmap-title 'markers' $data $out/heatmap.pdf

#$pegasus plot qc_violin --qc-type gene --cluster-labels $group --attribute Donor --subplot-size 7,5 --qc-xtick-font 5 --qc-line-width 0.5 $data  $out/qc_violin.pdf


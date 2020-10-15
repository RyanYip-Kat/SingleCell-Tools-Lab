R=/home/ye/anaconda3/envs/r-base/bin/Rscript
r1=/home/ye/Work/R/Cytof/dr2020/src/Make_Counts.R
r2=/home/ye/Work/R/Cytof/dr2020/src/Create_SingleCellExperiment.R
r3=/home/ye/Work/R/Cytof/dr2020/src/sceToseurat.R
r4=/home/ye/Work/R/Cytof/dr2020/src/antigen_plot.R

path=$1
number=$2
out=$3

echo "### Step 1"
$R $r1 --batch_correct --path $path --number $number --outdir $out

echo "### Step 2"
$R $r2 --count "$out/counts.rds" --outdir $out

echo "### Step 3"
$R $r3 --hm "$out/harmony.rds" --sce "$out/sce.rds" --outdir $out

echo "### Setp 4"
$R $r4 --seurat "$out/seurat.rds" --outdir "$out/antigen_plot" 

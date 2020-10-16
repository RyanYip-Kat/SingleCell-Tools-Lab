
#/home/ye/anaconda3/envs/scanpy/bin/pyscenic ctx  loom/aging-11/BC/adj.tsv \
#       	/home/ye/Work/R/SCENIC/pyScenic/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather  \
#	--output loom/aging-11/BC/reg.csv \
#	--num_workers 12 \
#	--annotations_fname "/home/ye/Work/R/SCENIC/pyScenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
#        --expression_mtx_fname loom/aging-11/BC/scenic.loom \
#        --mode custom_multiprocessing \
#        --mask_dropouts	

adj=$1
out=$2
loom=$3
/home/ye/anaconda3/envs/scanpy/bin/pyscenic ctx $adj \
        /home/ye/Work/R/SCENIC/pyScenic/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather  /home/ye/Work/R/SCENIC/pyScenic/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather \
        --output $out/reg.csv \
        --num_workers 12 \
        --annotations_fname "/home/ye/Work/R/SCENIC/pyScenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
        --expression_mtx_fname $loom \
        --mode custom_multiprocessing \
        --mask_dropouts


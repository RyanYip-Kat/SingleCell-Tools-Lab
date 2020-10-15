docker run -itd -v "/home/ye/Work/Python/SingleCell/Project/Isoform:/Isoform" \
    gregoryschwartz/too-many-cells:2.0.0.0 make-tree \
    --matrix-path /Isoform/NGS_union_matrix/k12_14_15/total/matrix.csv \
    --labels-file /Isoform/NGS_union_matrix/k12_14_15/total/label.csv \
    --draw-collection "PieRing" \
    --output /Isoform/NGS_union_matrix/k12_14_15/total/out \
    > NGS_clusters.csv


#docker run -it --rm -v "/home/ye/Work/Python/SingleCell/Project/Isoform:/Isoform" \
#    gregoryschwartz/too-many-cells:2.0.0.0 make-tree \
#    --prior /Isoform/NGS_union_matrix/not_cd14_vim/out \
#    --labels-file /Isoform/NGS_union_matrix/not_cd14_vim/label.csv \
#    --smart-cutoff 4 \
#    --min-size 1 \
#    --draw-collection "PieRing" \
#    --output /Isoform/NGS_union_matrix/not_cd14_vim/out_pruned \
#    > NGS_clusters_pruned.csv

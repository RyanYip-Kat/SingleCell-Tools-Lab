docker run -itd -v "/home/ye/Work/Python/SingleCell/Project/Isoform:/Isoform" \
    gregoryschwartz/too-many-cells:2.0.0.0 make-tree \
    --matrix-path /Isoform/TGS_union_matrix/min_cells_3/G1/6h/matrix \
    --labels-file /Isoform/TGS_union_matrix/min_cells_3/G1/6h/label.csv \
    --draw-collection "PieChart" \
    --draw-node-number \
    --output /Isoform/TGS_union_matrix/min_cells_3/G1/6h/out \
    > TGS_clusters.csv

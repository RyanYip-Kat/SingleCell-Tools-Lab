docker run -it --rm  -v "/home/ye/Work/Python/SingleCell/Project/Isoform:/Isoform" gregoryschwartz/too-many-cells:2.0.0.0 differential --matrix-path /Isoform/TGS_union_matrix/min_cells_3/G1/3h/matrix --prior /Isoform/TGS_union_matrix/min_cells_3/G1/3h/out --nodes "([4,6,7],[11,12,13])" > differential_3h_k14_k15.csv


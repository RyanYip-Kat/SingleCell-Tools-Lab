import csv
import os
import scipy.io
 
# peak-bc matrix
 
matrix_dir = "/opt/sample345/outs/filtered_peak_bc_matrix"
mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx"))
 
peaks_path = os.path.join(human_matrix_dir, "peaks.bed")
peaks = [(row[0], int(row[1]), int(row[2])) for row in csv.reader(open(peaks_path), delimiter="\t")]
 
barcodes_path = os.path.join(human_matrix_dir, "barcodes.tsv")
barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]

 
# tf-bc matrix
 
matrix_dir = "/opt/sample345/outs/filtered_tf_bc_matrix"
mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx"))
 
motifs_path = os.path.join(human_matrix_dir, "motifs.tsv")
motif_ids = [row[0] for row in csv.reader(open(motifs_path), delimiter="\t")]
motif_names = [row[1] for row in csv.reader(open(motifs_path), delimiter="\t")]


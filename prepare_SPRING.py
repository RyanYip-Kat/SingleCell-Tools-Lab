import numpy as np
import scipy

import argparse
import os

import sys
sys.path.append("./SPRING_dev/data_prep")

from spring_helper import *
from doublet_detector import *
from collections import defaultdict

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--path",type=str,default=None)

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

input_path=args.path
if os.path.isfile(input_path + '/matrix.npz'):
    E = scipy.sparse.load_npz(input_path + '/matrix.npz')
else:
    E = scipy.io.mmread(input_path + '/matrix.mtx').T.tocsc()
    scipy.sparse.save_npz(input_path + '/matrix.npz', E, compressed=True)

print E.shape

gene_list = np.array(load_genes(input_path + 'genes.tsv', delimiter='\t', column=1))
total_counts = E.sum(1).A.squeeze()
E = tot_counts_norm(E)[0]


main_spring_dir = args.outdir
np.savetxt(main_spring_dir + '/genes.txt', gene_list, fmt='%s')
np.savetxt(main_spring_dir + '/total_counts.txt', total_counts)


print 'Saving hdf5 file for fast gene loading...'
save_hdf5_genes(E, gene_list, main_spring_dir + '/counts_norm_sparse_genes.hdf5')

##############
print 'Saving hdf5 file for fast cell loading...'
save_hdf5_cells(E, main_spring_dir + '/counts_norm_sparse_cells.hdf5')

##############
save_sparse_npz(E, main_spring_dir + '/counts_norm.npz', compressed = False)


t0 = time.time()

save_path = main_spring_dir + '/full'

out = make_spring_subplot(E, gene_list, save_path,
                    normalize = False, tot_counts_final = total_counts,
                    min_counts = 3, min_cells = 3, min_vscore_pctl = 85,show_vscore_plot = True,
                    num_pc = 30,
                    k_neigh = 4,
                    num_force_iter = 500)

np.save(save_path + '/cell_filter.npy', np.arange(E.shape[0]))
np.savetxt(save_path + '/cell_filter.txt',  np.arange(E.shape[0]), fmt='%i')

print 'Finished in %i seconds' %(time.time() - t0)

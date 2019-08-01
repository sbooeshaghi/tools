#!/usr/bin/env python3
import sys
sys.path.insert(0, r'/home/sina/tools/scripts')

import matrix_utils as mu
import plotting_utils as pu
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pickle
import os
import datetime

def make_comparison(a_folder, b_folder, c_folder):
    A = mu.import_adata(a_folder)
    B = mu.import_adata(b_folder, cr=True)
    C = mu.import_adata(c_folder, cr=True)

    # Sort lexigraphically
    A.obs_names_make_unique()
    # A.var_names_make_unique()

    B.obs_names_make_unique()
    # B.var_names_make_unique()

    A = A[:, A.var.sort_index().index]
    A = A[A.obs.sort_index().index]

    B = B[B.obs.sort_index().index]
    B = B[:, B.var.sort_index().index]


    A = mu.basic_process(A)
    B = mu.basic_process(B)
    C = mu.basic_process(C)


    joint, common, A_var, B_var = mu.barcode_sets(A, B)

    A, B = mu.filter_adata(A, B, obs=True, var=True)

    Af, C = mu.filter_adata(A, C, by_C=True)
    Bf, C = mu.filter_adata(B, C, by_C=True)

    A_AB, M_AB = mu.MA(Af, Bf)

    cc_raw = mu.T_sparse_M_corr(A[:,common.index].layers['log1p'].T, B[:,common.index].layers['log1p'].T)
    cc_filtered = mu.T_sparse_M_corr(Af.layers['log1p'].T, Bf.layers['log1p'].T)

    Af = mu.compute_tsvd(Af)
    Bf = mu.compute_tsvd(Bf)

    Af = mu.compute_tsne(Af)
    Bf = mu.compute_tsne(Bf)

#    dist_AA, dist_AB = mu.l1_dist(Af.layers['log1p'].T, Bf.layers['log1p'].T)

    return A, B, Af, Bf, cc_raw, cc_filtered, A_AB, M_AB, joint, common, A_var, B_var

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make the comparsion between kallisto and Cell Ranger.')
    parser.add_argument("-a", type=str, help="kallisto folder")
    parser.add_argument("-b", type=str, help="Cell Ranger raw folder")
    parser.add_argument("-c", type=str, help="Cell Ranger filtered folder")
    parser.add_argument("outdir", type=str, help="Output for saving files.")

    args = parser.parse_args()

    A, B, Af, Bf, cc_raw, cc_filtered, A_AB, M_AB, joint, common, A_var, B_var = make_comparison(args.a, args.b, args.c)


    dataset_shortname = args.a.split("/")[-2]
    savedir = os.path.join(args.outdir + dataset_shortname + "/")
    print("[{}] Saving files..".format(datetime.datetime.now()))
    try:
        os.mkdir(save_dir)
    except:
        print("\t Folder already made.")


    A.write(savedir + 'A.h5ad')
    B.write(savedir + 'B.h5ad')
    Af.write(savedir + 'Af.h5ad')
    Bf.write(savedir + 'Bf.h5ad')

    with open(savedir + 'A_var.pkl', 'wb') as handle:
        pickle.dump(A_var, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(savedir + 'B_var.pkl', 'wb') as handle:
        pickle.dump(B_var, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(savedir + 'common.pkl', 'wb') as handle:
        pickle.dump(common, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(savedir + 'joint.pkl', 'wb') as handle:
        pickle.dump(joint, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # with open(savedir + 'dist_AA.pkl', 'wb') as handle:
    #     pickle.dump(dist_AA, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # with open(savedir + 'dist_AB.pkl', 'wb') as handle:
    #     pickle.dump(dist_AB, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(savedir + 'cc_raw.pkl', 'wb') as handle:
        pickle.dump(cc_raw, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(savedir + 'cc_filtered.pkl', 'wb') as handle:
        pickle.dump(cc_filtered, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(savedir + 'M_AB.pkl', 'wb') as handle:
        pickle.dump(M_AB, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(savedir + 'A_AB.pkl', 'wb') as handle:
        pickle.dump(A_AB, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print("[{}] Files have been saved to {}".format(datetime.datetime.now(), savedir))


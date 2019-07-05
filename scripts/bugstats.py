#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np


def umicollision(bug):
    print("Making unique bug file..")
    print("Done.")

    print("Grouping the genes by bcs, umi..")
    gb = ubug.groupby(["bcs", "umi"])[["gene"]].nunique()
    print("Done.")
    # distribution of number of nonunique genes per unique umi (in one cell) summed across all cells.
    collision = gb["gene"].values

    print("Getting a list of the bad umis..")
    dupes = gb[gb["gene"]>1]
    dupes.index = dupes.index.droplevel(0)

    # umi sequence and count: if the umi collided between genes more than once in a cell 
    # sum the number of collisions across all cells for that umi. The resulant number is the count.
    badumis = dupes.groupby("umi")[["gene"]].sum()


    return collision, badumis

def multiplicities(bug):
    v = bug["mul"].value_counts()
    return v



if __name__=="__main__":
    parser = argparse.ArgumentParser(description="processes bug file and saves (1) distribution of between gene collisions (2) list of umnis that collide between genes more than once (3) value counts of mul to feed into preseq.")
    parser.add_argument("--b", help="bug file in txt format")

    parser.add_argument("outdir", type=str, help="output directory for bug file")
    args = parser.parse_args()

    print("Loading files..")

    bug = pd.read_csv(args.b, header=None, names=["bcs" ,"umi", "gene", "mul"], sep="\t")
    print("Making unique bug..")
    ubug = bug[bug["gene"].map(lambda l: "," not in str(l))]
    print("Done")
    
    print("File loaded, computing umi collisions..")
    collision, badumis = umicollision(bug)
    print("Done.")

    print("Getting multiplicity value counts")
    sat = multiplicities(bug)
    print("Done.")

    badumis.to_csv(args.outdir + "badumis.txt", index=True, header=False, sep="\t")
    np.savetxt(args.outdir + "collision.txt", collision, fmt='%i')
    np.savetxt(args.outdir + "value_counts_gene.txt", np.c_[pd.Series(sat).sort_index().index,pd.Series(sat).sort_index().values], fmt='%i')
    
    with open(args.outdir + "k.txt", "w") as f:
        f.write("{}".format(ubug.shape[0]))

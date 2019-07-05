#!/usr/bin/env python

import sys
import argparse
import pandas as pd
import numpy as np
from scipy.optimize import fsolve


def prob(N, k):
    return (1-(1-1/N)**k-k/N*(1-1/N)**(k-1))
def solve(v, *data):
    Ne, k = v
    c, ng, num_unique_umis = data
    s = 0
    for i in ng:
        s += i*prob(Ne, i)

    return (Ne*prob(Ne, k) - s - int(k*c), num_unique_umis - Ne*(1-(1-1/Ne)**(k)))

def estN(bug):
    Nest=[]
    kest=[]
    ubug = bug[bug["gene"].map(lambda l: "," not in str(l))]
    umipercell = ubug.groupby("bcs")[["umi"]].count()
    umipercellpergene = ubug.groupby(["bcs", "gene"])[["umi"]].count()
    counter = 0
    for cellwewant in umipercell.umi.nlargest(100).keys()[0:50]:
        cellbug = ubug[ubug["bcs"] == cellwewant]
        num_unique_umis = cellbug.umi.nunique()
        gb = cellbug.groupby(["bcs", "umi"])[["gene"]].nunique()
        collision = gb["gene"].values

        pofbwgenecol = collision[collision>1].shape[0]/collision.shape[0]

        ng = umipercellpergene.loc[cellwewant].umi.values
        expcol = pofbwgenecol

        Ntmp, ktmp = fsolve(solve, (1000000, 10000), args=(expcol, ng, num_unique_umis))
        Nest.append(Ntmp)
        kest.append(ktmp)
        counter += 1
#        ",".join(map(str, Nest))
    return [int(i) for i in Nest], [int(i) for i in kest]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="estimate N, number of umis from which we are sampling.")
    parser.add_argument("--b", help="bug file in text format")
    parser.add_argument("outdir", type=str, help="output directory for Nest.txt")
    args = parser.parse_args()

    print("Loading bug file..")
    bug = pd.read_csv(args.b, header=None, names=["bcs", "umi", "gene", "mul"], sep="\t")

    Nest, kest = estN(bug)
    
    Nest = ",".join(map(str, Nest))
    kest = ",".join(map(str, kest))

    with open(args.outdir + "Nest.txt", "w") as f:
        f.write("{}".format(Nest))


    with open(args.outdir + "kest.txt", "w") as f:
        f.write("{}".format(kest))

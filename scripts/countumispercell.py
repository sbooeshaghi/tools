#!/usr/bin/env python

import argparse
import pandas as pd

def countumis(bug):
   gb = bug.groupby("bcs")["umi"].nunique()
   return ",".join(map(str, gb))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Count the number of unique umis per cell. Returns list..")
    parser.add_argument("--b", help="bug file in text format")
    parser.add_argument("outdir", type=str, help="output directory for umispercell.txt")
    args = parser.parse_args()
    
    print("Loading bug file..")
    bug = pd.read_csv(args.b, header=None, names=["bcs", "umi", "gene", "mul"], sep="\t")
    print("Done.")
    
    print("Counting UMIs....")
    umispercell = countumis(bug)
    print("Saving..")
    with open(args.outdir + "umispercell.txt", "w") as f:
        f.write("{}".format(umispercell))

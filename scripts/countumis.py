#!/usr/bin/env python

import argparse
import pandas as pd

def countumis(bug):
   gb = bug["umi"].nunique()
   return gb 

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Count the number of unique umis.")
    parser.add_argument("--b", help="bug file in text format")
    parser.add_argument("outdir", type=str, help="output directory for Nest.txt")
    args = parser.parse_args()
    
    print("Loading bug file..")
    bug = pd.read_csv(args.b, header=None, names=["bcs", "umi", "gene", "mul"], sep="\t")
    print("Done.")
    
    print("Counting UMIs....")
    N = countumis(bug)
    print("Saving..")
    with open(args.outdir + "N.txt", "w") as f:
        f.write("{}".format(N))

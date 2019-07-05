#!/usr/bin/env python

import argparse
import pandas as pd

def getsuppI(bug):
    ubug = bug[bug["gene"].map(lambda l: "," not in str(l))]
    gb = ubug.groupby(["bcs", "gene"])[["umi"]].nunique() 
    return gb 

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="get the umi counts per cell per gene to estimate counts lost.")
    parser.add_argument("--b", help="bug file in text format")
    parser.add_argument("outdir", type=str, help="output directory for Nest.txt")
    args = parser.parse_args()
    
    print("Loading bug file..")
    bug = pd.read_csv(args.b, header=None, names=["bcs", "umi", "gene", "mul"], sep="\t")
    print("Done.")
    
    suppI = getsuppI(bug)
    print("Saving..")
    suppI.to_csv(args.outdir + "suppI.txt", header=None, sep="\t")

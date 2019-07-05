#!/usr/bin/env python
import os
from scipy import sparse, io
from glob import glob
import scanpy as sc
fs = glob("/home/sina/projects/bus/velocity/bus_out/brain_1M_out/brain_1M_barcodes/kallisto*/")
def main():
    count = 0
    print("Writing file order..")
    with open("fileorder.txt", "w") as f:
        for i in fs:
            f.write(i)
            f.write("\n")

    print("Loading initial matrices..")
    s = io.mmread(fs[0] + "spliced.mtx")
    u = io.mmread(fs[0] + "unspliced.mtx")
    count += 1
    print("Stacking matrices..")
    for i in fs[1:]:
        count+=1
        spliced=io.mmread(i + "spliced.mtx")
        unspliced=io.mmread(i + "unspliced.mtx")
        s = sparse.vstack([s, spliced])
        u = sparse.vstack([u, unspliced])
        print("Finished {} of 133".format(count))
    print("Done")
    return s, u

if __name__ == "__main__":
    sX, uX = main()
    s = sc.AnnData(sX)
    u = sc.AnnData(uX)
    s.write("spliced.h5")
    u.write("unspliced.h5")

#for f in 
#spliced = io.mmread()
#unspliced = io.mmread()
#s = sparse.vstack([s, spliced])
#u = sparse.vstack([u, unspliced])

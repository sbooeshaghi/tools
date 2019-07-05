#!/usr/bin/env python
import pandas as pd
from scipy import io, sparse
import numpy as np

fnames = "/home/sina/projects/bus/velocity/bus_out/brain_1M_out/fileorder.txt"
cnt = 0
with open(fnames, "r") as f:
    for line in f.readlines()[::-1]:
        cnt += 1
        
        savedir = line[0:-1]
        line = ("/".join(line.split("/")[0:7]) + "/brain_1M/" + "/".join(line.split("/")[9::]))[0:-1]
        crf = "/".join(line.split("/")[0:8]) + "/cellranger_" + "_".join(line.split("/")[8].split("_")[1::])
        # print("k line: ", line)
        # print("savedir:", savedir)
        # print("c line: " , crf)
        # print("\n") 
        # print("Loading barcodes..")
        sbcs = pd.read_csv(line + "/spliced.barcodes.txt", header=None)
        ubcs = pd.read_csv(line + "/unspliced.barcodes.txt", header=None)
        cbcs = pd.read_csv(crf + "/barcodes.tsv" , header=None)
        cbcs[0] = cbcs[0].str.slice(0, 16) 

        print("CR Shape: ", cbcs.shape)
        print("U  Shape: ", sbcs.shape)
        print("S  Shape: ", ubcs.shape)

        print("Filtering barcodes by cellranger..")
        sidx = sbcs[0].isin(cbcs[0]).values
        uidx = ubcs[0].isin(cbcs[0]).values
       
        sidx2 = np.where(sidx)[0]
        uidx2 = np.where(uidx)[0]

        print("Reading in sparse matrices..")
        sX = io.mmread(line + "/spliced.mtx").tocsr()
        uX = io.mmread(line + "/unspliced.mtx").tocsr()

        print("Filtering sparse matrices by filtered barcodes..")
        fsX = sX[sidx2,:]
        fuX = uX[uidx2,:]

        print("CR Shape: ", cbcs.shape)
        print("S  Shape: ", fsX.shape)
        print("U  Shape: ", fuX.shape)

        print("Saving filtered barcodes..")
        sbcs[sidx].to_csv(savedir + "/spliced.barcodes.txt", header=None, index=False)
        ubcs[uidx].to_csv(savedir + "/unspliced.barcodes.txt", header=None, index=False)

        print("Saving filtered sparse matrices..")
        io.mmwrite(savedir + "/spliced.mtx", fsX)
        io.mmwrite(savedir + "/unspliced.mtx", fuX)
        print("Finished {} out of 133..\n".format(cnt))

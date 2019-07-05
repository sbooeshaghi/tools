#!/usr/bin/env python

from operator import itemgetter
import sys
import json

def parse_bed(bed):
    with open(bed, "r") as f:
        # want to retain intron information
        # ie what transcript the intron belongs to
        d = {}
        for line in f:
            line = line.split("\t")
            chrom = line[0]
            if chrom not in d.keys():
                d[chrom] = []
            low = int(line[1])
            high = int(line[2])
            trn = line[3]
            trn = trn.split("_")[2]
            d[chrom].append((low, high, trn))
        return d
        
def merge_intervals(intervals):
    trn_list = []
    sorted_intervals = sorted(intervals, key=itemgetter(0))
    if not sorted_intervals:  # no intervals to merge
        return
    low, high, trn = sorted_intervals[0]
    trn_list.append(trn)
    for iv in sorted_intervals[1:]:
        if iv[0] <= high:
            high = max(high, iv[1])
            trn_list.append(iv[2])
        else:
            yield low, high, trn_list
            trn_list = []
            low, high , trn = iv
            trn_list.append(iv[2])
    yield low, high, trn_list

if __name__ == "__main__":
    bed_file = sys.argv[1] 
    d = parse_bed(bed_file)
    #print(json.dumps(d, indent = 4))
    with open("intervals.bed", "w") as f:
        for i in d.keys():
            for entry in merge_intervals(d[i]):
                f.write(i + "\t" + str(entry[0]) + "\t" + str(entry[1]) + "\t" + ",".join(entry[2]) + "\n")

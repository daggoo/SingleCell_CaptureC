import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

import numpy as np
import pandas
import h5py

import cooler
import os

folder="hg19/"

filenames = sorted([j for j in os.listdir(folder) if j.endswith(".cool") ])

for fl in filenames:
    print fl.split(".cool")[0]
    c = cooler.Cooler(folder + fl)
    #print c.info
    
    
    bins = c.bins()[:]
    counts={}
    
    for i in range (0, 100):
        counts[i]=0
    
    pix = c.pixels()  # select some pixels with unannotated bins
    
    #print len(pix)
    
    pix = list(pix[:]['count'])
    for i, entry in enumerate(pix):
        
        if entry in counts:
            counts[entry]+=1
        else:
            print entry
    
    
    for i in  range (0, 100):
        if counts[i]!=0:
            print str(i) + "s =" ,
            print str(counts[i]) + ";" ,

    print "sum = " ,
    print sum(counts.values())
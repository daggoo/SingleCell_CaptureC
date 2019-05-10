
print "1"
from hiclib.fragmentHiC import HiCdataset
print "2"
from mirnylib.systemutils import fmap,setExceptionHook
print "3"
from mirnylib.genome import Genome 

print "4.0"
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

print "4"
from mirnylib import plotting
print "5"
import numpy as np 
print "6"
import os
print "7"
import sys
print "8"
from defineGenome import getGenome
from mirnylib import h5dict
from hiclib import binnedData

setExceptionHook()
def ensure(f):
    d = os.path.dirname(f)
    if os.path.isdir(d):
        return f
    else:
        try:
            os.makedirs(d)
        except:
            raise ValueError("Cannot create directory")
    return f


coolerResolutions = [1000]




def refineDataset(filenames, create=True, delete=False, parseInMemory=True):
    """
    Parameters
    ----------

    filenames[0] is a list of filenames of incoming files
    filenames[1] is a folder for outgoing file
    filenames[2] is a working genome, that is output directory
    filenames[3] is an enzyme for a given experiment


    create : bool, optional
        If True, parse each file.
        If False, assume that files were already parsed
        (e.g. if you are just playing around with filtering parameters)
    delete : bool, optional
        If True, delete parsed files after merging.
        Man, these files may be huge... if you don't have a 10TB RAID, this may be useful.
    parseInMemory : bool, optional
        Perform parsing input files in memory.

    """
    in_files = filenames[0]
    out_file = filenames[1]

    statFolder = os.path.join("statistics", out_file)

    workingGenome = filenames[2]
    enzyme = filenames[3]

    if create == True:  # if we need to parse the input files (.hdf5 from mapping).
        def parse_onename(onename):
            np.random.seed()
            #Parsing individual files
            if parseInMemory == True:
                finalname = onename + "_parsed.frag"
                #if not os.path.exists(finalname):
                if True:

                    #create dataset in memory, parse and then save to destination
                    TR = HiCdataset("bla" + str(np.random.randint(100000000000)), genome=getGenome(workingGenome),
                                    maximumMoleculeLength=500,enzymeName = enzyme,tmpFolder = "tmp",
                                    inMemory=True)  # remove inMemory if you don't have enough RAM

                    TR.parseInputData(dictLike=onename)
                    folder = os.path.split(onename)[0]
                    print(onename)
                    TR.save(ensure(finalname))
                    folder, fname = os.path.split(onename)
                    statSubFolder = os.path.join(statFolder, folder)
                   

                    TR.printMetadata(saveTo=ensure(os.path.join(statSubFolder, fname + ".stat")))
                else:
                    print("skipping parsed: ", onename)
            else:
                #Create dataset at destination, parse on HDD, then no need to save.
                TR = HiCdataset(ensure(onename + "_parsed.frag"),
                                genome=getGenome(workingGenome),enzymeName = enzyme,tmpFolder = "tmp",
                                maximumMoleculeLength=500, mode='w')
                TR.parseInputData(dictLike=onename, enzymeToFillRsites=enzyme)
                TR.printMetadata(saveTo=ensure(os.path.join(statFolder, onename + ".stat")))
        list(map(parse_onename, in_files))
        "Merging files alltogether, applying filters"
        TR = HiCdataset(ensure(out_file + "_merged.frag"),
                        genome=getGenome(workingGenome),enzymeName = enzyme,tmpFolder = "tmp",dictToStoreIDs="h5dict",
                        mode="w")
        TR.merge([i + "_parsed.frag" for i in in_files])
            #Merge in all parsed files from one experiment

        if delete == True:  # cleaning up parsed files
            for delFile in [i + "_parsed.frag" for i in in_files]:
                os.remove(delFile)

        "Now opening new dataset for refined data, and performing all the filtering "
        TR = HiCdataset(out_file + "_refined.frag",enzymeName = enzyme,
                        genome=getGenome(workingGenome),tmpFolder = "tmp",dictToStoreIDs="h5dict",
                        mode='w')
        TR.load(out_file + "_merged.frag")

        #----------------------------Set of filters applied -------------
        TR.filterDuplicates()
        #TR.save(out_file+".dat")
        #TR.filterExtreme(cutH=0.0001, cutL=0)
        #TR.filterRsiteStart()
        #TR.filterLarge()
        TR.writeFilteringStats()
        TR.printMetadata(saveTo=statFolder + ".stat")

        #------------------------End set of filters applied----------

    else:
        #If merging & filters has already been done, just load files
        TR = HiCdataset(out_file + "_working.frag",enzymeName = enzyme,
                        mode='w', genome=getGenome(workingGenome))
        TR.load(out_file + "_refined.frag")
        TR.printMetadata(saveTo=statFolder + ".stat")

    print("----->Building Raw heatmap at different resolutions")
    TR.printStats()
    for res in coolerResolutions: 
        TR.saveCooler(out_file+".{0}.cool".format(res), res)
        #pass
        #TR.saveHeatmap(out_file+".{0}.hm".format(res), res)
         


#This code is actually parsing datasets.tsv file

try:
    filename = sys.argv[1]
except:
    filename = "datasets_Golov.tsv"

fsplit = os.path.split(filename)
if len(fsplit[0]) > 0:
    os.chdir(fsplit[0])
filename = fsplit[1]

lines = open(filename).readlines()
lines = [i for i in lines if i[0] != "#"]
lines = [i.split() for i in lines if (len(i) > 3) and (i[0] != "#")]
for line in lines:
    if len(line) != 5:
        print("bad line", line)
        raise

for i in lines:
    if len(i) == 4:
        print(i) 
        #Now we assume that enzyme is fifth argument in datasets.tsv or second commandline argument
        try:
            i.append(sys.argv[2])
        except:
            print("bad line", i)
            print("Please specify enzyme as a second command line argument, or as a fifth column")
            raise ValueError("Please specify enzyme as a second command line argument, or as a fifth column")

assert False not in [len(i) == 5 for i in lines]

dataFiles = lines
#experiment is defined by experiment name, replica name, genome and enzyme 
experimentNames = set((i[1], i[2], i[3], i[4]) for i in dataFiles)
print "\n\n======================print experimentNames=======================\n\n"

print experimentNames
print "\n\n"

byExperiment = []
combinedExperimentNames = []

for experiment in experimentNames:
    workingGenome = experiment[2]
    enzyme = experiment[3]
    filenames = [i[0] for i in dataFiles if (i[1], i[2], i[3], i[4]) == experiment]
    outName = "{0}-{1}-{3}".format(*experiment)

    #Filenames, path to save, genome, enzyme
    byExperiment.append((filenames, os.path.join(workingGenome, outName), workingGenome, enzyme))
    


#Now running refineDataset for each experiment
for i in byExperiment:
    refineDataset(i, create=True, delete=False, parseInMemory=False)
    pass

'''
raw_heatmap = h5dict.h5dict('hg19/SS-R1-HindIII.1000000.hm', mode='r') 
resolution = int(raw_heatmap['resolution'])

genomeName="hg19"
genome_db = getGenome(genomeName)

BD = binnedData.binnedData(resolution, genome_db)

BD.simpleLoad('hg19/SS-R1-HindIII.1000000.hm', 'HindIII')

plotting.plot_matrix(np.log(BD.dataDict['HindIII']))
plt.savefig("hg19/SS-R1-HindIII.1000000.png", dpi=300, figsize=(16, 16))
'''

import atexit
import glob
import os
import sys
import gzip

import logging
from hiclib import mapping
from mirnylib import h5dict, genome
import pickle
import numpy as np
logging.basicConfig(level=logging.DEBUG)
from mirnylib.systemutils import setExceptionHook
from defineGenome import getGenome

def cleanDirectory(dirName):
    for i in os.listdir(dirName):
        os.remove(os.path.join(dirName, i))
        
def cleanFile(filename):
    if os.path.exists(filename):
        os.remove(filename)

def calculateStep(length, minlen, approxStep=10, maxSteps=4):
    """returns minimum length and step based on the
    length of sequence and proposed minimum length"""

    actualDif = length - minlen
    if actualDif < approxStep * 0.6:
        return length, 100

    numIter = np.array(np.around(actualDif / float(approxStep)), dtype=int)
    if numIter == 0:
        numIter = 1
    if numIter > maxSteps:
        numIter = maxSteps
    actualStep = actualDif / numIter

    minlen = length - actualStep * numIter

    return minlen, actualStep

def splitSingleFastq(filename, outFile, splitBy=4000000, convertReadID=lambda x:x):

    inFile = os.path.abspath(filename)

    pread = subprocess.Popen(["gunzip", inFile, "-c"],
                             stdout=subprocess.PIPE, bufsize=-1)
    inStream = pread.stdout

    halted = False
    counters = []
    for counter in range(100000):

        outProc1 = gzipWriter(outFile.format(counter))
        outStream1 = outProc1.stdin

        for j in range(splitBy):

            line = inStream.readline()

            try:
                assert six.indexbytes(line,0) == 64 #"@"
            except AssertionError:
                print('Not fastq')
                print("bad line: {0}".format(line))                
                raise IOError("File is not fastq: {0}".format(filename))
            except IndexError:
                halted = True
                counters.append(j)
                break

            fastq_entry = (convertReadID(line), inStream.readline(),
                           inStream.readline(), inStream.readline())
            outStream1.writelines(fastq_entry)

        outProc1.communicate()
        print("finished block number", counter)

        if halted:
            if (counters[-1] < splitBy / 3) and (len(counters) > 1):
                f1 = outFile.format(counter - 1)
                f2 = outFile.format(counter)
                os.system("cat {0} {1} > {0}_tmp".format(f1, f2))
                shutil.move(f1 + "_tmp", f1)
                os.remove(f2)
                last = counters.pop()
                counters[-1] = counters[-1] + last
            print("Read counts", counters)
            return counters
        counters.append(splitBy)



def doOne(inData):
    file1, file2, outfile = inData
    print("Mapping {0} and {1} into {2}".format(*inData))

    
    for onefile in file1, file2:
        '''
        a = gzip.open(onefile, 'r')
        a.readline()
        length = len(a.readline()) - 1
        if length < 10:
            raise ValueError("Length of your sequence is {0}. Something is wrong".format(length))
        seqSkipStart = 0
        
        minlen, step = calculateStep(length - seqSkipStart, minMapLen)
        
        
        print "minlen" ,
        print minlen
        #minlen=15
        print minlen
        print step
        '''
        
        minlen = 15
        step = 10
        seqSkipStart = 0
        
        mapping.iterative_mapping(
            bowtie_path=bowtiePath,
            bowtie_index_path=bowtieIndex,
            fastq_path=onefile,
            out_sam_path=os.path.join(samFolder, os.path.split(onefile)[1] + ".sam"),
            seq_start=seqSkipStart,
            min_seq_len=minlen,  # for bacteria mimimal mappable length is 15 bp, so I start with something slightly longer
            len_step=step,  # and go with a usualy step
            nthreads=threads,  # on intel corei7 CPUs 4 threads are as fast as
                         # 8, but leave some room for you other applications
            # max_reads_per_chunk = 10000000,  #optional, on low-memory machines
            temp_dir=tmpDir,
            bowtie_flags=bowtieFlags,
            drop_sequences=False
            )

    os.remove(file1)
    os.remove(file2)
    
    
    # Second step. Parse the mapped sequences into a Python data structure,
    #    assign the ultra-sonic fragments to restriction fragments.
    mapped_reads = h5dict.h5dict(outfile)
    sf1, sf2 = [os.path.join(samFolder, os.path.split(onefile)[1] + ".sam") for onefile in [file1, file2]]
    
    print "\n\n=========sf========\n\n"
    print sf1
    print sf2
    print "\n\n"
    print "\n\n=========outfile========\n\n"
    print outfile
    print "\n\n"        
    

    mapping.parse_sam(sam_basename1=sf1, sam_basename2=sf2,
        out_dict=mapped_reads, genome_db=genome_db, save_seqs=False, maxReads=int(chunkSize*1.6), IDLen=50)
    for i in os.listdir(samFolder):
        if (os.path.split(file1)[1] in i) or (os.path.split(file2)[1] in i):
            #print("not deleting", i)
                        
            print("deleting", i)
            os.remove(os.path.join(samFolder, i))
            



#
#
os.system("rm -rf mapped-hg19")
os.system("rm -rf tmp")
os.system("mkdir -m 777 tmp")
        
mode="fastq"

inFastqDir="raw_fastq/"
tmpDir = "tmp/"
chunkSize = 10000000
sidePrefixes = ("_R1","_R2")  # a prefix preceeding .fastq.gz, which will be used to distinguish side 1 and side 2
genomeName = "hg19"
#genomeName = "/home/ubuntu/data3/HiCdata/DM_DpnII_Oct18/dm3"
threads = 8
bowtieIndex = "/home/ubuntu/rDNA_data/index/hg19/hg19mannual"
#bowtiePath = "/usr/bin/bowtie2" 
bowtiePath = "/home/ubuntu/tools/miniconda2/envs/main/bin/bowtie2" 
bowtieFlags = "--very-sensitive"

#bowtieFlags = "-D 20 -R 3 -N 0 -L 7 -i S,1,0.10"

bowtieFlags = "--very-sensitive --n-ceil L,0,0.1"

seqSkipStart = 0  # skip first 2 bp of the read, if you want
#minMapLen = 25  # start mapping at this length
genome_db = getGenome(genomeName)

if not os.path.exists(bowtiePath): raise

if mode == "sra":
    iterList = 2 * GEOids
elif mode == "fastq":
    iterList =  sorted(os.listdir(inFastqDir))

print iterList

for i in  iterList:
    if mode == "sra":
        sraNum = i
        expName = "SRR{0}".format(i)
        i = expName
        num = int(i[3:])
    if mode == "fastq":
        if not i.endswith(sidePrefixes[0] + ".fastq.gz"):
            if not i.endswith(sidePrefixes[1] + ".fastq.gz"):
                print("ignoring file", i, "does not end with fastq.gz")
            continue
        expName = i.replace(sidePrefixes[0] + ".fastq.gz", "")
        
    fastqFolder = os.path.join(tmpDir, "{0}-fastq-{1}".format(expName, genomeName))
    samFolder = os.path.join(tmpDir, "{0}-sams-{1}".format(expName, genomeName))
    
    print expName
    print fastqFolder
    print samFolder    
    
    savePath = "mapped-{0}".format(genomeName)
    saveFolder = os.path.join(savePath, expName)

    for folder in [samFolder, saveFolder, fastqFolder]:
        if not os.path.exists(folder):
            os.makedirs(folder)
            

    lockName = saveFolder + ".lock"
    completedName = os.path.join(saveFolder, "completed")
    
    if os.path.exists(completedName) and not os.path.exists(lockName):
        print("skipping", expName)
        continue
    
    print lockName
    if os.path.exists(lockName):
        print("someone is working on", expName)
        continue
    if os.path.exists(completedName) and os.path.exists(lockName):
        raise

    
    cleanDirectory(fastqFolder)
    cleanDirectory(samFolder)
    atexit.register(cleanDirectory, fastqFolder)
    atexit.register(cleanDirectory, samFolder)
    

    lock = open(lockName, "w")
    lock.close()
    atexit.register(cleanFile, lockName)
    
    
    firstSide = os.path.join(inFastqDir, expName + sidePrefixes[0] + ".fastq.gz")
    secondSide = os.path.join(inFastqDir, expName + sidePrefixes[1] + ".fastq.gz")
    
    
    counter = mapping.splitSingleFastq(firstSide, os.path.join(fastqFolder, expName + "_chunk{0:04d}_fhtagn_side1.fastq.gz"), chunkSize)
    mapping.splitSingleFastq(secondSide, os.path.join(fastqFolder, expName + "_chunk{0:04d}_fhtagn_side2.fastq.gz"), chunkSize)
    
    inFiles = [os.path.join(fastqFolder, i) for i in os.listdir(fastqFolder)]
    print inFiles
    
    for i in inFiles:
        atexit.register(cleanFile, i)

    inFiles1 = sorted([i for i in inFiles if "fhtagn_side1" in i])
    inFiles2 = sorted([i for i in inFiles if "fhtagn_side2" in i])
    assert len(inFiles1) == len(inFiles2)
    outFiles = [os.path.join(saveFolder, "chunk{0:04d}.hdf5".format(i + 1)) for i in range(len(inFiles1))]
    
    print outFiles   

    
    if len(inFiles1) == 0:
        raise ValueError("No files supplied")
    print list(zip(inFiles1, inFiles2, outFiles))
    list(map(doOne, list(zip(inFiles1, inFiles2, outFiles))))
    a = open(completedName, 'w')
    a.close()
    
    os.remove(lockName)    
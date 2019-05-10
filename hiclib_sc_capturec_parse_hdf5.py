import h5py    
import numpy as np    
from defineGenome import getGenome
import mirnylib
from mirnylib import h5dict, genome
from hiclib.fragmentHiC import HiCdataset
import cooler

from mirnylib.systemutils import fmap,setExceptionHook

from mirnylib import h5dict
from hiclib import binnedData

import os


import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

print "4"
from mirnylib import plotting



enzyme = "HaeIII"

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


path = "/home/ubuntu/RNAdata/Golov/hic/mapped-hg19/35_AGAGGCCT_reads_merged_len_adj/"
f1 = h5py.File(path + "chunk0001.hdf5" ,'r+')  

f1 = mirnylib.h5dict.h5dict(path + "chunk0001.hdf5" ,'r+')  


print f1["misc"]["genome"]["idx2label"]
chrm_conversion_table = f1["misc"]["genome"]["idx2label"]


genome_db = getGenome("hg19")
workingGenome = "hg19"

for key in genome_db.idx2label:
    print key ,
    print genome_db.idx2label[key] ,

print("Keys: %s" % f1.keys())

chrms1_key = list(f1.keys())[0]
chrms2_key = list(f1.keys())[1]

cuts1_key = list(f1.keys())[2]
cuts2_key = list(f1.keys())[3]

misc_key = list(f1.keys())[4]

strand1_key = list(f1.keys())[5]
strand2_key = list(f1.keys())[6]

print cuts1_key

'''

cuts1 = list(f1[cuts1_key])
cuts2 = list(f1[cuts2_key])

chrms1 = list(f1[chrms1_key])
chrms2 = list(f1[chrms2_key])

strand1 = list(f1[strand1_key])
strand2 = list(f1[strand2_key])

print len (cuts1)

print cuts1[0:9]
print cuts2[0:9]

print strand1[0:9]
print strand2[0:9]

print chrms1[0:9]
print chrms2[0:9]


count_valid = 0
count_duplicates = 0

table1 = {}
table2 = {}
table = {} 
with open(path + "dict_in_text", "w") as outfile:
    #outfile.write(misc)
    for i in range(len(cuts1)):
        if cuts1[i]!=-1 and cuts2[i]!=-1:
            count_valid+=1
        else:    
            id1 = str(cuts1[i]) + "-" + str(cuts2[i])
            chr1 = str(chrms1[i]) + "-" + str(chrms2[i])
            
            id2 = str(cuts1[i]) + "-" + str(cuts2[i])
            chr2 = str(chrms1[i]) + "-" + str(chrms2[i])
            
            if id1 in table:
                if chr1 in table[id1]:
                    count_duplicates+=1
                else:
                    table[id1].append(chr1)
            else:
                table[id1]=[chr1]
                
            if id2 in table:
                if chr2 in table[id2]:
                    count_duplicates+=1
                else:
                    table[id2].append(chr2)
            else:
                table[id2]=[chr2]     
            
        outfile.write(str(chrms1[i]) + "\t")
        outfile.write(str(chrms2[i]) + "\t")
        outfile.write(str(cuts1[i]) + "\t")
        outfile.write(str(cuts2[i]) + "\t")
        outfile.write(str(strand1[i]) + "\t")
        outfile.write(str(strand2[i]) + "\n") 
        
print count_duplicates 

'''
'''

in_files = ["mapped-hg19/35_AGAGGCCT_reads_merged_len_adj/chunk0001.hdf5"]

outName = "35_test"


workingGenome = "hg19"

out_file = os.path.join(workingGenome, outName)

statFolder = os.path.join("statistics", out_file)


enzyme = "HaeIII"

onename = in_files[0]
finalname = onename + "_parsed.frag"


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

print out_file + "_merged.frag"
"Merging files alltogether, applying filters"
TR = HiCdataset(ensure(out_file + "_merged.frag"),
                genome=getGenome(workingGenome),enzymeName = enzyme,tmpFolder = "tmp",dictToStoreIDs="h5dict",
                mode="w")
TR.merge([i + "_parsed.frag" for i in in_files])
    #Merge in all parsed files from one experiment


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


filename = onename + "_filtered.frag"
TR.save(ensure(finalname))

'''



'''
path = "/home/ubuntu/RNAdata/Golov/hic/hg19/"

f1 = mirnylib.h5dict.h5dict(path + "35_test_refined.frag" ,'r+')  


#print f1["misc"]["genome"]["idx2label"]

print("Keys: %s" % f1.keys())

chrms1_key = list(f1.keys())[0]
chrms2_key = list(f1.keys())[1]

cuts1_key = list(f1.keys())[2]
cuts2_key = list(f1.keys())[3]

misc_key = list(f1.keys())[4]

strand1_key = list(f1.keys())[5]
strand2_key = list(f1.keys())[6]


cuts1 = list(f1[cuts1_key])
cuts2 = list(f1[cuts2_key])

chrms1 = list(f1[chrms1_key])
chrms2 = list(f1[chrms2_key])

strand1 = list(f1[strand1_key])
strand2 = list(f1[strand2_key])

print len (cuts1)

print cuts1[0:9]
print cuts2[0:9]

print strand1[0:9]
print strand2[0:9]

print chrms1[0:9]
print chrms2[0:9]


count_valid = 0
count_duplicates1 = 0
count_duplicates2 = 0
count_order = 0

table1 = {}
table2 = {}
table = [] 
with open(path + "dict_in_text", "w") as outfile:
    #outfile.write(misc)
    for i in range(len(cuts1)):
        if cuts1[i] > cuts2[i]:
            count_order+=1
            
            
            
            id1 = str(cuts1[i]) + "-" + str(cuts2[i]) + "-" + str(chrms1[i]) + "-" + str(chrms2[i])            
            
            id2 = str(cuts2[i]) + "-" + str(cuts1[i]) + "-" + str(chrms2[i]) + "-" + str(chrms1[i])
            
            
            if id1 in table:                
                count_duplicates1+=1                
            else:
                table.append(id1)
                
            if id2 in table:
                
                    count_duplicates2+=1
    
            
        outfile.write("chr" + str(chrm_conversion_table[chrms1[i]]) + "\t")
        outfile.write("chr" + str(chrm_conversion_table[chrms2[i]]) + "\t")
        outfile.write(str(cuts1[i]) + "\t")
        outfile.write(str(cuts2[i]) + "\t")
        outfile.write(str(strand1[i]) + "\t")
        outfile.write(str(strand2[i]) + "\n") 
        
print count_order
print count_duplicates1
print count_duplicates2
'''


'''
path = "/home/ubuntu/RNAdata/Golov/hic/1attempt/hg19/"

#filenames = [j for j in os.listdir(path) if j.endswith("-R1-HaeIII_refined.frag")]
#filenames = [j for j in os.listdir(path) if j.endswith("SC35-R1-HaeIII_merged.frag")]
filenames = [j for j in os.listdir(path) if j.endswith("_refined.frag")]


for fl in sorted(filenames):
    
    
    print fl
    
    f1 = HiCdataset(path + fl.split("_refined.frag")[0] + "_refined_wo_man_dupl.frag",enzymeName = enzyme,
                    genome=getGenome(workingGenome),tmpFolder = "tmp",dictToStoreIDs="h5dict",
                    mode='w')
    f1.load(path + fl)
    
    
    #f1 = mirnylib.h5dict.h5dict(path + fl ,'r+')  
    
    
    print len(f1.chrms1)
    print len(f1.chrms2)
    print len(f1.cuts1)
    print len(f1.cuts2)
    
    chrms1_key = list(f1.keys())[0]
    chrms2_key = list(f1.keys())[1]
    
    cuts1_key = list(f1.keys())[2]
    cuts2_key = list(f1.keys())[3]
    
    misc_key = list(f1.keys())[4]
    
    strand1_key = list(f1.keys())[5]
    strand2_key = list(f1.keys())[6]
    
    
    cuts1 = list(f1[cuts1_key])
    cuts2 = list(f1[cuts2_key])
    
    chrms1 = list(f1[chrms1_key])
    chrms2 = list(f1[chrms2_key])
    
    strand1 = list(f1[strand1_key])
    strand2 = list(f1[strand2_key])
    
    
    
  
    
    count_valid = 0
    count_duplicates1 = 0
    count_duplicates2 = 0
    count_order = 0
    duplicates1 = []
    duplicates2 = []
    table = [] 
    
    table_coord = {}
    dup = 0
    count_wo_dup = 0
    def get_eps(x):        
        return range(x-5, x+5)
    
    eps = 5
    with open(path + fl + "_table", "w") as outfile ,\
         open(path + fl + "_table_wo_dupl", "w") as outfile_wo_dup:
        for i in range(len(cuts1)):
            dup = 0
            coord1 = cuts1[i]
            coord2 = cuts2[i]
            
            chr1 = str(chrms1[i])
            chr2 = str(chrms2[i]) 
            
            if cuts1[i] > cuts2[i]:
                count_order+=1                
                
            
                
            id1 = str(coord1) + "-" + str(coord2) + "-" + str(chr1) + "-" + str(chr2)   
            id2 = str(coord2) + "-" + str(coord1) + "-" + str(chr2) + "-" + str(chr1)
            
            
            
            #if id1 in table:
            #    count_duplicates1+=1 
            #else:
            #    table.append(id1)
            
            #if id2 in table:                    
            #    count_duplicates2+=1
            
            
            
            if chr1 in table_coord:
                for k in range(coord1 - eps, coord1 + eps + 1):
                    if k in table_coord[chr1] and chr2 in table_coord[chr1][k]:
                        for j in range(coord2 - eps, coord2 + eps + 1):
                            
                            id_j1 = str(k) + "-" + str(j) + "-" + str(chr1) + "-" + str(chr2) 
                                                        
                            if j in table_coord[chr1][k][chr2]:
                                count_duplicates1+=1 
                                dup = 1
                                if id_j1 not in duplicates1 and id_j1 not in duplicates2: 
                                    duplicates1.append(id_j1)                                
            
                            
            if chr2 in table_coord:
                for j in range(coord2 - eps, coord2 + eps + 1):
                    if j in table_coord[chr2] and chr1 in table_coord[chr2][j]:
                        for k in range(coord1 - eps, coord1 + eps + 1):   
                            
                            id_j2 = str(j) + "-" + str(k) + "-" + str(chr2) + "-" + str(chr1) 
                            
                            if chr2 in table_coord \
                               and j in table_coord[chr2] \
                               and chr1 in table_coord[chr2][j] \
                               and k in table_coord[chr2][j][chr1] :
                                
                                count_duplicates2+=1             
                                dup = 1                                  
                                                             
                                if id_j2 not in duplicates1 and id_j2 not in duplicates2:
                                    duplicates2.append(id_j2)                                 
                                               
                                    
                                

                
            

                                
            if chr1 not in table_coord:
                table_coord[chr1] = {}
            if coord1 not in table_coord[chr1]:
                table_coord[chr1][coord1] = {}
            if chr2 not in table_coord[chr1][coord1]:
                table_coord[chr1][coord1][chr2] = {}
            table_coord[chr1][coord1][chr2][coord2] = 1
                                
            outfile.write("chr" + str(chrm_conversion_table[chrms1[i]]) + "\t")
            outfile.write("chr" + str(chrm_conversion_table[chrms2[i]]) + "\t")
            outfile.write(str(cuts1[i]) + "\t")
            outfile.write(str(cuts2[i]) + "\t")
            outfile.write(str(strand1[i]) + "\t")
            outfile.write(str(strand2[i]) + "\n")        
            if dup == 0:
                outfile_wo_dup.write("chr" + str(chrm_conversion_table[chrms1[i]]) + "\t")
                outfile_wo_dup.write("chr" + str(chrm_conversion_table[chrms2[i]]) + "\t")
                outfile_wo_dup.write(str(cuts1[i]) + "\t")
                outfile_wo_dup.write(str(cuts2[i]) + "\t")
                outfile_wo_dup.write(str(strand1[i]) + "\t")
                outfile_wo_dup.write(str(strand2[i]) + "\n") 
                count_wo_dup+=1
            for k in sorted(duplicates1):
                pass
                #print k
                
    print fl ,
    print len (cuts1) ,
    print count_wo_dup ,
    #print count_order
    print len(duplicates1) ,
    print count_duplicates1 ,
    print len(duplicates2) ,
    print count_duplicates2    

'''



path = "/home/ubuntu/RNAdata/Golov/hic/1attempt/hg19/"

#filenames = [j for j in os.listdir(path) if j.endswith("-R1-HaeIII_refined.frag")]
#filenames = [j for j in os.listdir(path) if j.endswith("SC35-R1-HaeIII_merged.frag")]
filenames = [j for j in os.listdir(path) if j.endswith("_refined.frag")]

genomeName="hg19"
genome_db = getGenome(genomeName)


for fl in sorted(filenames):
    
    
    print fl
    
    out_file = path + fl.split("_refined.frag")[0] + "_refined_wo_man_dupl"
    
    f1 = HiCdataset(path + fl.split("_refined.frag")[0] + "_refined_wo_man_dupl.frag",enzymeName = enzyme,
                    genome=getGenome(workingGenome),tmpFolder = "tmp",dictToStoreIDs="h5dict",
                    mode='w')
    f1.load(path + fl)
    
    
    #f1 = mirnylib.h5dict.h5dict(path + fl ,'r+')  
    
    
    print len(f1.chrms1)
    
  
    
    cuts1 = list(f1.cuts1)
    cuts2 = list(f1.cuts2)
    
    chrms1 = list(f1.chrms1)
    chrms2 = list(f1.chrms2)
    
    strand1 = list(f1.strands1)
    strand2 = list(f1.strands2)
    
    
    cuts1_new = []
    cuts2_new = []
    
    chrms1_new = []
    chrms2_new = []
    
    strand1_new = []
    strand2_new = []
    

        
    
    
    
    '''
    print len(f1.chrms1)
    print len(f1.chrms2)
    print len(f1.cuts1)
    print len(f1.cuts2)    
    
        
    for name in f1.vectors:
        print "n1 " ,
        print name 
    '''
    dupl = np.zeros((f1.N, 2), dtype="int64", order="C")
    
    
    
    count_valid = 0
    count_duplicates1 = 0
    count_duplicates2 = 0
    count_order = 0
    duplicates1 = []
    duplicates2 = []
    table = [] 
    
    table_coord = {}
    dup = 0
    count_wo_dup = 0
    def get_eps(x):        
        return range(x-5, x+5)
    
    eps = 5
    
    
    with open(path + fl + "_table", "w") as outfile ,\
         open(path + fl + "_table_wo_dupl", "w") as outfile_wo_dup:
       
        for i in range(len(cuts1)):
            flag_duplicates = 0
            dup = 0
            coord1 = cuts1[i]
            coord2 = cuts2[i]
            
            chr1 = str(chrms1[i])
            chr2 = str(chrms2[i]) 
            
            if cuts1[i] > cuts2[i]:
                count_order+=1                
                
            
                
            id1 = str(coord1) + "-" + str(coord2) + "-" + str(chr1) + "-" + str(chr2)   
            id2 = str(coord2) + "-" + str(coord1) + "-" + str(chr2) + "-" + str(chr1)
            

           
            if chr1 in table_coord:
                for k in range(coord1 - eps, coord1 + eps + 1):
                    if k in table_coord[chr1] and chr2 in table_coord[chr1][k]:
                        for j in range(coord2 - eps, coord2 + eps + 1):
                            
                            id_j1 = str(k) + "-" + str(j) + "-" + str(chr1) + "-" + str(chr2) 
                                                        
                            if j in table_coord[chr1][k][chr2]:
                                count_duplicates1+=1 
                                dup = 1
                                if id_j1 not in duplicates1 and id_j1 not in duplicates2: 
                                    duplicates1.append(id_j1)                                
                                    
                            
            if chr2 in table_coord:
                for j in range(coord2 - eps, coord2 + eps + 1):
                    if j in table_coord[chr2] and chr1 in table_coord[chr2][j]:
                        for k in range(coord1 - eps, coord1 + eps + 1):   
                            
                            id_j2 = str(j) + "-" + str(k) + "-" + str(chr2) + "-" + str(chr1) 
                            
                            if chr2 in table_coord \
                               and j in table_coord[chr2] \
                               and chr1 in table_coord[chr2][j] \
                               and k in table_coord[chr2][j][chr1] :
                                
                                count_duplicates2+=1             
                                dup = 1                                  
                                                             
                                if id_j2 not in duplicates1 and id_j2 not in duplicates2:
                                    duplicates2.append(id_j2)                                 
                                               
                                    
                                

                
            

                                
            if chr1 not in table_coord:
                table_coord[chr1] = {}
            if coord1 not in table_coord[chr1]:
                table_coord[chr1][coord1] = {}
            if chr2 not in table_coord[chr1][coord1]:
                table_coord[chr1][coord1][chr2] = {}
                
            table_coord[chr1][coord1][chr2][coord2] = 1
                                
            outfile.write("chr" + str(chrm_conversion_table[chrms1[i]]) + "\t")
            outfile.write("chr" + str(chrm_conversion_table[chrms2[i]]) + "\t")
            outfile.write(str(cuts1[i]) + "\t")
            outfile.write(str(cuts2[i]) + "\t")
            outfile.write(str(strand1[i]) + "\t")
            outfile.write(str(strand2[i]) + "\n")        
            if dup == 0:
                
                outfile_wo_dup.write("chr" + str(chrm_conversion_table[chrms1[i]]) + "\t")
                outfile_wo_dup.write("chr" + str(chrm_conversion_table[chrms2[i]]) + "\t")
                outfile_wo_dup.write(str(cuts1[i]) + "\t")
                outfile_wo_dup.write(str(cuts2[i]) + "\t")
                outfile_wo_dup.write(str(strand1[i]) + "\t")
                outfile_wo_dup.write(str(strand2[i]) + "\n") 
                
                
                count_wo_dup+=1
                
                cuts1_new.append(cuts1[i])
                cuts2_new.append(cuts2[i])
                
                chrms1_new.append(chrms1[i])
                chrms2_new.append(chrms2[i])
                
                strand1_new.append(strand1[i])
                strand2_new.append(strand2[i])

                
    print fl ,
    print len (cuts1) ,
    print count_wo_dup ,
    #print count_order
    print len(duplicates1) ,
    print count_duplicates1 ,
    print len(duplicates2) ,
    print count_duplicates2    
    
    f1._setData("cuts1", cuts1_new)
    f1._setData("cuts2", cuts2_new)
    f1._setData("chrms1", chrms1_new)
    f1._setData("chrms2", chrms2_new)
    f1._setData("strands1", strand1_new)
    f1._setData("strands2", strand2_new)
    
    f1.N = len (strand2_new)
    
    f1.filterDuplicates()
    f1.writeFilteringStats()
    f1.printMetadata(saveTo=out_file + ".stat")
    
    for res in [1000]: 
        f1.saveCooler(out_file+".{0}.cool".format(res), res)
        fl_cooler = out_file+".{0}.cool".format(res)    
        c = cooler.Cooler(fl_cooler)
        
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
        
        
    for res in [1000000]: 
        #TR.saveCooler(out_file+".{0}.cool".format(res), res)
        #pass
        f1.saveHeatmap(out_file+".{0}.hm".format(res), res)   

        
        BD = binnedData.binnedData(res, genome_db)
        
        BD.simpleLoad(out_file+".{0}.hm".format(res), enzyme)
        
        plotting.plot_matrix(BD.dataDict[enzyme], clip_min = 0, clip_max = 1,  cmap='Blues')
        
        
        #plt.colorbar(extend='both')
        #plt.clim(0, 1);
        #, label = "'viridis'"
        plt.savefig(out_file+".{0}.png".format(res), dpi=300, figsize=(16, 16))    
        
        
        
        plt.close()
        
        for res in [1000000, 200000]: 
            f1.saveCooler(out_file+".{0}.cool".format(res), res)
            f1.saveHeatmap(out_file+".{0}.hm".format(res), res)  
        
        
        
          

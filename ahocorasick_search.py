import ahocorasick
import cPickle
import os
import gzip
import gc
import sys
import fileinput
import re 
import itertools
import operator 


path = ""
barcodesFile = "barcodes.txt"
barcodes_by_idx = {}
barcodes_by_seq = {}


complDict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', "N" : "N"}

def returnReverseCompliment(s):
    revcomp = ""
    for nt in s:
        revcomp = complDict[nt] + revcomp
    return revcomp
#
#
def returnCompliment(s):
    comp = ""
    for nt in s:
        comp = comp + complDict[nt]
    return comp
#
#
def returnReverse(s):
    rev = ""
    for nt in s:
        rev = nt + rev 
    return rev
#
#
def create_automaton():
    barcodes_string_fw = ""
    barcodes_string_rv = ""
    with open(barcodesFile, 'r') as bc:
        for i, l in enumerate(bc):
            #l = "GGA" + l.rstrip() + "GATCGA-TCC"
            l_fw = "" + l.rstrip() + "GATCGA"
            l_rv = returnReverseCompliment(l_fw)
            
            rcl=returnReverseCompliment(l.rstrip())
            #barcodes_string_fw = barcodes_string_fw + l  + " " + rcl + " "   
            barcodes_string_fw = barcodes_string_fw + l_fw  + " "   
            barcodes_string_rv = barcodes_string_rv + l_rv  + " "   
            
    automaton_fw = ahocorasick.Automaton()   
    automaton_rv = ahocorasick.Automaton()  
    
    with open(barcodesFile.split(".txt")[0] + "_indexed", 'w+') as barcodes_indexed:
        
        for idx, key in enumerate(barcodes_string_fw.rstrip().split(" ")):
            automaton_fw.add_word(key, (idx * 2, key))
            barcodes_indexed.write(str(idx * 2) + "\t" + key + "\n")  
        for idx, key in enumerate(barcodes_string_rv.rstrip().split(" ")):
            automaton_rv.add_word(key, (idx * 2 + 1, key))
            barcodes_indexed.write(str(idx * 2 + 1) + "\t" + key + "\n")   
            
    automaton_fw.make_automaton()
    automaton_rv.make_automaton()
    with open(barcodesFile.split(".txt")[0] + "_automaton_fw_pickled", 'w+') as barcodes_automaton_file:
        cPickle.dump(automaton_fw, barcodes_automaton_file)
    with open(barcodesFile.split(".txt")[0] + "_automaton_rv_pickled", 'w+') as barcodes_automaton_file:
        cPickle.dump(automaton_rv, barcodes_automaton_file)    
#
#
def load_automaton_from_pickle():
    with open(barcodesFile.split(".txt")[0] + "_automaton_fw_pickled", 'r') as bc:
        B_fw = cPickle.load(bc)  
    with open(barcodesFile.split(".txt")[0] + "_automaton_rv_pickled", 'r') as bc:
        B_rv = cPickle.load(bc)        
    return (B_fw, B_rv)    

#
#
def sorting_with_automaton(folderName, B_fw, B_rv):
    dirname=folderName.split(".")[0]
    if not os.path.isdir(dirname):
        os.system("mkdir " + dirname)    
    resultingFiles={}
    
    for el in B_fw:        
        idx = B_fw.get(el)[0]
        if idx % 2 == 0:
            resultingFiles[str(idx) + "_" + el[0:8] + "_single"] = open (dirname + "/" + str(idx) + "_" + el[0:8] + "_single", 'w+')    
        if idx % 2 == 1:            
            resultingFiles[str(idx) + "_" + el[6:] + "_single"] = open (dirname + "/" + str(idx) + "_" + el[6:] + "_single", 'w+') 
    
    
    
    
    
    for el in B_rv:        
        idx = B_rv.get(el)[0]
        if idx % 2 == 0:
            resultingFiles[str(idx) + "_" + el[0:8] + "_single"] = open (dirname + "/" + str(idx) + "_" + el[0:8] + "_single", 'w+')    
        if idx % 2 == 1:            
            resultingFiles[str(idx) + "_" + el[6:] + "_single"] = open (dirname + "/" + str(idx) + "_" + el[6:] + "_single", 'w+') 

    print sorted(resultingFiles.keys())
    
    count=0  
    iterator_fw=B_fw.iter("")
    iterator_rv=B_rv.iter("")
    count_to_report=0
    for lines in itertools.izip_longest(*[fileinput.input()]*4):
        found=False
        count+=4
        if (count % 10000000) == 0: 
            print count            
        line=lines[1].rstrip()
        iterator_fw.set(line, True)
        
        for end_index, (idx, original_value) in iterator_fw:
            found = True
            if idx % 2 == 0:
                fl_name = str(idx) + "_" + original_value[0:8] + "_single"                          
            start_index = end_index - len(original_value) + 1
            #assert line[start_index:start_index + len(original_value)] == original_value
            to_write=str(((start_index - 3, end_index + 3, (idx, original_value)))) + "\t" + line 
            resultingFiles[fl_name].write(to_write)   
            resultingFiles[fl_name].write("\t" + lines[0]) 
            resultingFiles[fl_name].write(lines[1]) 
            resultingFiles[fl_name].write(lines[2]) 
            resultingFiles[fl_name].write(lines[3]) 
            if count_to_report == 0:
                break
        if found == False:
            iterator_rv.set(line, True)
            for end_index, (idx, original_value) in iterator_rv:                
                if idx % 2 == 1:            
                    fl_name = str(idx) + "_" + original_value[6:] + "_single"             
                start_index = end_index - len(original_value) + 1
                #assert line[start_index:start_index + len(original_value)] == original_value
                to_write=str(((start_index - 3, end_index + 3, (idx, original_value)))) + "\t" + line 
                resultingFiles[fl_name].write(to_write)   
                resultingFiles[fl_name].write("\t" + lines[0]) 
                resultingFiles[fl_name].write(lines[1]) 
                resultingFiles[fl_name].write(lines[2]) 
                resultingFiles[fl_name].write(lines[3])
                if count_to_report == 0:
                    break                
                
    for fl in resultingFiles.values():
        fl.close()
#    
#
def download_from_s3():
    folders=["results/strict2nt/Ulianov-1-13_R1_001", "results/strict2nt/Ulianov-1-13_R2_001"]    
    with open("barcodes_indexed", 'r') as bcs:
        for folder in folders:
            for _, l in enumerate(bcs):
                l=l.rstrip().split("\t")
                command="wget -P " \
                    + path + folder + " -O " + path + folder + "/" + l[0] + "_" + l[1] + "_single" \
                    + " https://s3.us-east-2.amazonaws.com/magnumdata/Golov/strict2nt/" \
                    + folder.split("/")[2] + "/" + l[0] + "_" + l[1] + "_single"
                os.system(command)    
    
#
#
def cut_read(read_data, adapters, foldername):
    results_1st=[]
    results_2nd=[]
    m={}
    count1=0
    count2=0
    count3=0
    total_adapter_number=0
    
    read=read_data[1].rstrip()
    
    pair_id = 0
    
    #print read_data[0].rstrip().split("@")
    #print read_data[0].rstrip().split("@")[1].split(" ")    
    #print "@" + read_data[0].rstrip().split("@")[1].split(" ")[0]
    
    
    name1 = "@" + read_data[0].rstrip().split("@")[1].split(" ")[0]
    name2 = read_data[0].rstrip().split("@")[1].split(" ")[1][1:]
    
    quality=read_data[3].rstrip()
    for adp in adapters.keys():
        m[adp]=[mt.start() for mt in re.finditer(adapters[adp], read) ]
    
    if 'fwD' in adapters.keys():
        #print 'fwD'
        for st in m['fwD']:
            
            total_adapter_number+=1
            #results.append(read[0:st])
            #results.append(read[st + len(adapters['fwD']) + 1 : ])
            r1 = returnReverse(read[0:st])
            r2 = read[st + len(adapters['fwD']) + 1 : ]
            
            if len(r1) > 25 and len(r2) > 25:
                count1+=1
                id1 = name1 + ":" + foldername + ":PiD" + str(pair_id) + " 1" + name2  
                qual1 = returnReverse(quality[0:st])
                results_1st.append([id1, r1, "+", qual1])
                            
                id2 = name1 + ":" + foldername + ":PiD" + str(pair_id) + " 2" + name2  
                qual2 = returnReverse(quality[st + len(adapters['fwD']) + 1 : ])
                pair_id+=1
                results_2nd.append([id2, r2, "+", qual2])
                
            if len(r1) < 15 or len(r2) < 15:    
                count2+=1
            if len(r1) < 25 or len(r2) < 25:   
                count3+=1
        if 'fwS' in adapters.keys():
            #print 'fwS'
            for st in m['fwS']:
                if st not in m['fwD']:
                    total_adapter_number+=1
                    #results.append(read[0:st])
                    #results.append(read[st + len(adapters['fwS']) + 1 : ])   
                    r1=returnReverse(read[0:st])
                    r2=read[st + len(adapters['fwS']) + 1 : ]          
                    
                    if len(r1) > 15 and len(r2) > 15:
                        count1+=1
                        id1 = name1 + ":" + foldername + ":PiD" + str(pair_id) + " 1" + name2  
                        qual1 = returnReverse(quality[0:st])
                        results_1st.append([id1, r1, "+", qual1])
                                                                      
                        id2 = name1 + ":" + foldername + ":PiD" + str(pair_id) + " 2" + name2  
                        qual2 = returnReverse(quality[st + len(adapters['fwS']) + 1 : ])
                        pair_id+=1
                        results_2nd.append([id2, r2, "+", qual2])                        
                    if len(r1) < 15 or len(r2) < 15:    
                        count2+=1
                    if len(r1) < 25 or len(r2) < 25:   
                        count3+=1                    
        if 'rvS' in adapters.keys():
            #print 'rvS'
            for st in m['rvS']:
                if (st - len (adapters['rvS'])) not in m['fwD']:
                    total_adapter_number+=1
                    #results.append(read[0:st])
                    #results.append(read[st + len(adapters['rvS']) + 1 : ])     
                    r1=returnReverse(read[0:st])                    
                    r2=read[st + len(adapters['rvS']) + 1 : ]  
                    
                    if len(r1) > 15 and len(r2) > 15:
                        count1+=1                        
                        id1 = name1 + ":" + foldername + ":PiD" + str(pair_id) + " 1" + name2  
                        qual1 = returnReverse(quality[0:st])
                        results_1st.append([id1, r1, "+", qual1])
                                               
                        id2 = name1 + ":" + foldername + ":PiD" + str(pair_id) + " 2" + name2  
                        qual2 = returnReverse(quality[st + len(adapters['rvS']) + 1 : ])
                        pair_id+=1
                        results_2nd.append([id2, r2, "+", qual2])
                        
                    if len(r1) < 15 or len(r2) < 15:    
                        count2+=1
                    if len(r1) < 25 or len(r2) < 25:   
                        count3+=1                    
    else:
        if 'fwS' in adapters.keys():
            #print 'fwS'
            for st in m['fwS']:     
                total_adapter_number+=1
                #results.append(read[0:st])
                #results.append(read[st + len(adapters['fwS']) + 1 : ])        
                r1=returnReverse(read[0:st])   
                r2=read[st + len(adapters['fwS']) + 1 : ]
                
                if len(r1) > 15 and len(r2) > 15:
                    count1+=1
                    id1 = name1 + ":" + foldername + ":PiD" + str(pair_id) + " 1" + name2  
                    qual1 = returnReverse(quality[0:st])
                    results_1st.append([id1, r1, "+", qual1])
                                      
                    id2 = name1 + ":" + foldername + ":PiD" + str(pair_id) + " 2" + name2  
                    qual2 = returnReverse(quality[st + len(adapters['fwS']) + 1 : ])
                    pair_id+=1
                    results_2nd.append([id2, r2, "+", qual2])
                    
                if len(r1) < 15 or len(r2) < 15:    
                    count2+=1
                if len(r1) < 25 or len(r2) < 25:   
                    count3+=1                
        if 'rvS' in adapters.keys():
            #print 'rvS'
            for st in m['rvS']:   
                total_adapter_number+=1
                #results.append(read[0:st])
                #results.append(read[st + len(adapters['rvS']) + 1 : ])   
                r1=returnReverse(read[0:st])   
                r2=read[st + len(adapters['rvS']) + 1 : ]
                if len(r1) > 15 and len(r2) > 15:
                    count1+=1
                    id1 = name1 + ":" + foldername + ":PiD" + str(pair_id) + " 1" + name2  
                    qual1 = returnReverse(quality[0:st])
                    results_1st.append([id1, r1, "+", qual1])
                    
                    
                    id2 = name1 + ":" + foldername + ":PiD" + str(pair_id) + " 2" + name2  
                    qual2 = returnReverse(quality[st + len(adapters['rvS']) + 1 : ])
                    pair_id+=1
                    results_2nd.append([id2, r2, "+", qual2])
                    
                if len(r1) < 15 or len(r2) < 15:    
                    count2+=1
                if len(r1) < 25 or len(r2) < 25:   
                    count3+=1                
    return {'1st' : results_1st, '2nd' : results_2nd, 'count1' : count1, 'count2' : count2, 'count3' : count3, 'total_adapter_number' : total_adapter_number}
#
#
def parse_reads():
    barcodes_idx={}
    statistics={}
    #folders=["results/strict2nt/Ulianov-1-13_R1_001", "results/strict2nt/Ulianov-1-13_R2_001"]
    folders=["Ulianov-1-13_R1_001", "Ulianov-1-13_R2_001"]
    #folders=["Ulianov-1-13_R1_001"]
    #folders=["results/strict2nt/Ulianov-1-13_R2_001"]   
    with open("barcodes_indexed", 'r') as bcs:
        for i, l in enumerate(bcs):
            l = l.rstrip().split("\t")
            barcodes_idx[int(l[0])] = l[1]
    count1={}
    count2={}
    count3={}        
    count_total_adapter={}  
    for idx in barcodes_idx:  
        if idx % 2 != 0:
            continue      
        
            #pass
            #break
        statistics[idx]={}
        count1[idx]={}
        count2[idx]={}
        count3[idx]={}      
        count_total_adapter[idx]={} 
        for folder in folders:              
            if "_R1_" in folder:
                fsuffix="UR1"
            elif "_R2_" in folder:
                fsuffix="UR2"
            else:
                fsuffix="ERROR"
                
                
            if idx/2 not in [35, 37, 50, 58, 60, 85, 88, 89]:
                continue
    
                
            print idx/2 ,
            folderName = folder
            print folderName
            
            statistics[idx][folderName]={'fwS' : 0, 'rvS' : 0, 'fwD' : 0, 
                                         'fwSfwD' : 0, 'fwSrvS' : 0, 'fwDrvS' : 0,
                                         'fwSfwDrvS' : 0}
            fl_forward = path + folder + "/" + str(idx) + "_" + barcodes_idx[idx][0:8] + "_single" 
            fl_reverse = path + folder + "/" + str(idx + 1) + "_" + barcodes_idx[idx + 1][6:] + "_single" 
            
            fl_readsR1 = path + folder + "/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_readsR1" 
            fl_readsR2 = path + folder + "/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_readsR2" 
            
            fwS = "A" + barcodes_idx[idx] + "T"
            rvS = "A" + barcodes_idx[idx + 1] + "T"
            #fwD = barcodes_idx[idx][ : -2] + barcodes_idx[idx + 1][2 : ]
            fwD = "A" + barcodes_idx[idx] + "TA" + barcodes_idx[idx + 1] + "T"
            rvD = returnCompliment(fwD)
            
            count1[idx][folder]=0
            count2[idx][folder]=0
            count3[idx][folder]=0
            count_total_adapter[idx][folder]=0
            
            count_fwD=0
            count_fwS=0
            count_rvS=0
            #print fwD
            #print fwS
            
            reads_1st=[]
            reads_2nd=[]
            with open(fl_forward, 'r') as data:
                for lines in itertools.izip_longest(*[data]*4):
                    l=lines[1]
                    read=lines[1]
                    ln=l.rstrip().split("\t") 
                    #pos = re.findall('\d+', ln[0]) 
                    #key = re.findall('[A-Z]{20}', ln[0])[0] 
                    
                    count_fwD=l.count(fwD)
                    count_fwS=l.count(fwS)
                    count_rvS=l.count(rvS)
                    
                    #print count_fwD
                    #print count_fwS
                    #print count_rvS
                    if count_fwD == count_fwS and count_fwD == count_rvS and count_fwD > 0:
                        statistics[idx][folderName]['fwD'] = statistics[idx][folderName]['fwD'] + 1  
                        stat=cut_read(lines, {'fwD' : fwD}, fsuffix)
                        count1[idx][folder]+=stat['count1']
                        count2[idx][folder]+=stat['count2']
                        count3[idx][folder]+=stat['count3']
                        count_total_adapter[idx][folder]+=stat['total_adapter_number']
                        reads_1st.extend(stat['1st'])
                        reads_2nd.extend(stat['2nd'])
                        
                    elif count_fwD < count_fwS and count_fwD == count_rvS and count_fwD > 0:
                        statistics[idx][folderName]['fwSfwD'] = statistics[idx][folderName]['fwSfwD'] + 1
                        stat=cut_read(lines, {'fwS' : fwS,'fwD' : fwD}, fsuffix)
                        count1[idx][folder]+=stat['count1']
                        count2[idx][folder]+=stat['count2']
                        count3[idx][folder]+=stat['count3']  
                        count_total_adapter[idx][folder]+=stat['total_adapter_number']
                        reads_1st.extend(stat['1st'])
                        reads_2nd.extend(stat['2nd'])
                        
                    elif count_fwD == count_fwS and count_fwD < count_rvS and count_fwD > 0:
                        statistics[idx][folderName]['fwDrvS'] = statistics[idx][folderName]['fwDrvS'] + 1    
                        stat=cut_read(lines, {'fwD' : fwD, 'rvS' : rvS}, fsuffix)
                        count1[idx][folder]+=stat['count1']
                        count2[idx][folder]+=stat['count2']
                        count3[idx][folder]+=stat['count3']   
                        count_total_adapter[idx][folder]+=stat['total_adapter_number']
                        reads_1st.extend(stat['1st'])
                        reads_2nd.extend(stat['2nd'])
                        
                    elif count_fwD == 0  and count_fwS > 0 and count_rvS == 0:
                        statistics[idx][folderName]['fwS'] = statistics[idx][folderName]['fwS'] + 1     
                        stat=cut_read(lines, {'fwS' : fwS}, fsuffix)
                        count1[idx][folder]+=stat['count1']
                        count2[idx][folder]+=stat['count2']
                        count3[idx][folder]+=stat['count3']     
                        count_total_adapter[idx][folder]+=stat['total_adapter_number']
                        reads_1st.extend(stat['1st'])
                        reads_2nd.extend(stat['2nd'])
                        
                    elif count_fwD == 0  and count_fwS > 0 and count_rvS > 0:
                        statistics[idx][folderName]['fwSrvS'] = statistics[idx][folderName]['fwSrvS'] + 1  
                        stat=cut_read(lines, {'fwS' : fwS, 'rvS' : rvS}, fsuffix)
                        count1[idx][folder]+=stat['count1']
                        count2[idx][folder]+=stat['count2']
                        count3[idx][folder]+=stat['count3']    
                        count_total_adapter[idx][folder]+=stat['total_adapter_number']
                        reads_1st.extend(stat['1st'])
                        reads_2nd.extend(stat['2nd'])
                        
                    elif count_fwD > 0  and count_fwS > count_fwD and count_rvS > count_fwD:
                        statistics[idx][folderName]['fwSfwDrvS'] = statistics[idx][folderName]['fwSfwDrvS'] + 1  
                        stat=cut_read(lines, {'fwS' : fwS, 'fwD' : fwD, 'rvS' : rvS}, fsuffix)
                        count1[idx][folder]+=stat['count1']
                        count2[idx][folder]+=stat['count2']
                        count3[idx][folder]+=stat['count3']                        
                        count_total_adapter[idx][folder]+=stat['total_adapter_number']
                        reads_1st.extend(stat['1st'])
                        reads_2nd.extend(stat['2nd'])                        
                        
            with open(fl_reverse, 'r') as data:
                for lines in itertools.izip_longest(*[data]*4):
                    l=lines[1]
                    ln=l.rstrip().split("\t") 
                    read=lines[1]
                    #pos = re.findall('\d+', ln[0]) 
                    #key = re.findall('[A-Z]{20}', ln[0])[0] 
                    count_fwD=l.count(fwD)
                    count_fwS=l.count(fwS)
                    count_rvS=l.count(rvS)
                    
                    if count_fwS == 0 and count_fwD == 0 and count_rvS > 0:
                        statistics[idx][folderName]['rvS'] = statistics[idx][folderName]['rvS'] + 1  
                        stat=cut_read(lines, {'rvS' : rvS}, fsuffix)
                        count1[idx][folder]+=stat['count1']
                        count2[idx][folder]+=stat['count2']
                        count3[idx][folder]+=stat['count3']  
                        count_total_adapter[idx][folder]+=stat['total_adapter_number']
                        reads_1st.extend(stat['1st'])
                        reads_2nd.extend(stat['2nd'])                        
                    elif (count_fwS > 0 or count_fwD > 0) and count_rvS > 0:
                        print "Error1"
            with open (fl_readsR1, 'w+') as wr:
                for rd in reads_1st:
                    wr.write(rd[0] + "\n")
                    wr.write(rd[1] + "\n")
                    wr.write(rd[2] + "\n")
                    wr.write(rd[3] + "\n")
            with open (fl_readsR2, 'w+') as wr:
                for rd in reads_2nd:
                    wr.write(rd[0] + "\n")
                    wr.write(rd[1] + "\n")
                    wr.write(rd[2] + "\n")
                    wr.write(rd[3] + "\n")
                            
            '''
            print "count_fwS = " ,
            print statistics[idx][folderName]['fwS']
            print "count_fwD = " ,
            print statistics[idx][folderName]['fwD'] 
            print "count_rvS = " ,
            print statistics[idx][folderName]['rvS'] 
            print "count_fwSfwD = " ,
            print statistics[idx][folderName]['fwSfwD']
            print "count_fwDrvS = " ,
            print statistics[idx][folderName]['fwDrvS'] 
            print "count_fwSrvS = " ,
            print statistics[idx][folderName]['fwSrvS']               
            print "count_fwSfwDrvS = " ,
            print statistics[idx][folderName]['fwSfwDrvS']      
            
            print "count1 = " ,
            print count1 
            print "count2 = " ,
            print count2 
            print "count3 = " ,
            print count3 
            '''
    
    with open("stat2", 'w') as out:    
        out.write("idx" + "\t" + "barcode" + "\t" + "fwS" + "\t" + "fwD" + "\t" + "rvS" + "\t" + "fwSfwD"
                  + "\t" + "fwDrvS" + "\t" + "fwSrvS" + "\t" + "fwSfwDrvS" + "\t" + "count1" + "\t" + "count2" + "\t" + "count3" + "\t" + "count_total_number_of_adapters"  
                  + "\t" + "fwS" + "\t" + "fwD" + "\t" + "rvS" + "\t" + "fwSfwD"
                  + "\t" + "fwDrvS" + "\t" + "fwSrvS" + "\t" + "fwSfwDrvS" + "\t" + "count1" + "\t" + "count2" + "\t" + "count3" + "\t" + "count_total_number_of_adapters" 
                  + "\n")        
        for idx in sorted(statistics.keys()): 
            if idx/2 not in [35, 37, 50, 58, 60, 85, 88, 89]:
                continue
            
            out.write(str(idx/2) + "\t" )
            out.write("A" + barcodes_idx[idx] + "T" + "\t" )
            out.write(str(statistics[idx]["Ulianov-1-13_R1_001"]['fwS']) + "\t" )
            out.write(str(statistics[idx]["Ulianov-1-13_R1_001"]['fwD']) + "\t" )  
            out.write(str(statistics[idx]["Ulianov-1-13_R1_001"]['rvS']) + "\t" )
            out.write(str(statistics[idx]["Ulianov-1-13_R1_001"]['fwSfwD']) + "\t" )
            out.write(str(statistics[idx]["Ulianov-1-13_R1_001"]['fwDrvS']) + "\t" )  
            out.write(str(statistics[idx]["Ulianov-1-13_R1_001"]['fwSrvS']) + "\t" )
            out.write(str(statistics[idx]["Ulianov-1-13_R1_001"]['fwSfwDrvS']) + "\t" )
            out.write(str(count1[idx]["Ulianov-1-13_R1_001"]) + "\t" )
            out.write(str(count2[idx]["Ulianov-1-13_R1_001"]) + "\t" )
            out.write(str(count3[idx]["Ulianov-1-13_R1_001"]) + "\t" )
            out.write(str(count_total_adapter[idx]["Ulianov-1-13_R1_001"]) + "\t" )
            
            out.write(str(statistics[idx]["Ulianov-1-13_R2_001"]['fwS'])  + "\t" )
            out.write(str(statistics[idx]["Ulianov-1-13_R2_001"]['fwD'])  + "\t" )   
            out.write(str(statistics[idx]["Ulianov-1-13_R2_001"]['rvS'])  + "\t" )
            out.write(str(statistics[idx]["Ulianov-1-13_R2_001"]['fwSfwD']) + "\t" )
            out.write(str(statistics[idx]["Ulianov-1-13_R2_001"]['fwDrvS']) + "\t" )  
            out.write(str(statistics[idx]["Ulianov-1-13_R2_001"]['fwSrvS']) + "\t" )
            out.write(str(statistics[idx]["Ulianov-1-13_R2_001"]['fwSfwDrvS']) + "\t" )            
            out.write(str(count1[idx]["Ulianov-1-13_R2_001"]) + "\t" )
            out.write(str(count2[idx]["Ulianov-1-13_R2_001"]) + "\t" )
            out.write(str(count3[idx]["Ulianov-1-13_R2_001"]) + "\t" )            
            out.write(str(count_total_adapter[idx]["Ulianov-1-13_R2_001"]) + "\n" ) 



#
#
def concatenate_reads():    
    barcodes_idx={}
    
    folders=["Ulianov-1-13_R1_001", "Ulianov-1-13_R2_001"]    
    with open("barcodes_indexed", 'r') as bcs:
        for i, l in enumerate(bcs):
            l = l.rstrip().split("\t")
            barcodes_idx[int(l[0])] = l[1]
    count1={}
    count2={}
    count3={}        
    count_total_adapter={}  
    for idx in barcodes_idx:  
        if idx % 2 != 0:
            continue     
        
        if idx/2 not in [35, 37, 50, 58, 60, 85, 88, 89]:
            continue
        
        print idx/2 ,
        
        fl_readsR1_U1 = path + folders[0] + "/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_readsR1"  
        fl_readsR1_U2 = path + folders[1] + "/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_readsR1" 
        
        filenames = [fl_readsR1_U1, fl_readsR1_U2]                
        outputFile= path + "hic/raw_fastq/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_reads_merged_R1.fastq"        
        with open(outputFile, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)            
                        
        fl_readsR2_U1 = path + folders[0] + "/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_readsR2"  
        fl_readsR2_U2 = path + folders[1] + "/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_readsR2" 
        filenames = [fl_readsR2_U1, fl_readsR2_U2]                
        outputFile= path + "hic/raw_fastq/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_reads_merged_R2.fastq"        
        with open(outputFile, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)         

#
#
def check_duplicates():
    barcodes_idx={}
    
    folders=["Ulianov-1-13_R1_001", "Ulianov-1-13_R2_001"]    
    with open("barcodes_indexed", 'r') as bcs:
        for i, l in enumerate(bcs):
            l = l.rstrip().split("\t")
            barcodes_idx[int(l[0])] = l[1]
    count1={}
    count2={}
    count3={}        
    count_total_adapter={}  
    #for idx in barcodes_idx:
    
    prefix=[]
    
    for n1 in ["A", "T", "G", "C", "N"] :
        for n2 in ["A", "T", "G", "C", "N"]:
            for n3 in ["A", "T", "G", "C", "N"]:
                for n4 in ["A", "T", "G", "C", "N"]:
                    for n5 in ["A", "T", "G", "C", "N"]:
                        for n6 in ["A", "T", "G", "C", "N"]:
                            for n7 in ["A", "T", "G", "C", "N"]:
                                for n8 in ["A", "T", "G", "C", "N"]:
                                    for n9 in ["A", "T", "G", "C", "N"]:
                                        
                                            R1[n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8 + n9  ] ={}                    
                                            prefix.append(n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8 + n9 )
                                            
    for idx in [35 * 2, 50 * 2, 58 * 2, 60 * 2, 85 * 2, 88 * 2, 89 * 2]:
        if idx % 2 != 0:
            continue     
        
        if idx/2 not in [35, 37, 50, 58, 60, 85, 88, 89]:
            continue
        
        print "Idx = " ,
        print idx/2     
        
        fl_R1 = path + "hic/raw_fastq/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_reads_merged_R1.fastq"
        fl_R2 = path + "hic/raw_fastq/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_reads_merged_R2.fastq"
        
        fl_R1_rmdup = path + "hic/raw_fastq/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_reads_merged_rmdup_R1.fastq"
        fl_R2_rmdup = path + "hic/raw_fastq/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_reads_merged_rmdup_R2.fastq"
        
        fl_duplicates = path + "hic/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_ids_duplicates"
        
        R1 = {}
        R1_all = {}
        R2 = {}
        

                        
                
        duplicate_pairs={}
        duplicates_R1={}
        duplicates_R2={}
        duplicates=[]
        duplicates_by_prefix={}
        with open(fl_R1, 'r') as data:
            for lines in itertools.izip_longest(*[data]*4):
                id=lines[0].rstrip().split(" ")[0]
                read=lines[1].rstrip()                
                R1[read[0:9]][id]=read      
                R1_all[id]=read   
                
        with open(fl_R2, 'r') as data:
            for lines in itertools.izip_longest(*[data]*4):
                id=lines[0].rstrip().split(" ")[0]                
                read=lines[1].rstrip()               
                R2[id]=read
        
        
        count_r1=0
        count_r2=0
        count_missed=0
        
        count_total_rp=len(R2.keys())
        for id in R1:
            count_r1+=1
            if id not in R2:
                #print "Warning " + id
                count_missed+=1
                continue
                    
        #print count_r1         
        #print count_missed
        
        count=0
        count_dup=0
        for pr in prefix:
            
            
            sorted_reads = sorted(R1[pr].items(), key=operator.itemgetter(1))
            
            
            ids= [read[0] for read in sorted_reads]
            
            len_ids = len(ids)
            if len_ids == 0:
                continue
            print "prefix = " ,
            print pr ,
            
            print " number of reads = " ,
            print len_ids ,
            
            count +=  len_ids
            print " checked reads " ,
            print count ,
            print " out of " ,
            print count_total_rp 
            
            duplicates_local = []
            
            for i, id in enumerate(ids):
                if i % 1000 == 0:
                    pass
                    #print i
                
                
                #if id in duplicates:
                #    continue
                if id not in R2:
                    print "Warning " + id
                    break
                
                
                
                r1 = R1[pr][id]
                r2 = R2[id]
                
                
                for it in ids[i+1:]:
                    if it == id:
                        continue
                    #if it in duplicates_local:
                    #    continue                
                    if it not in R2:
                        print "Warning " + it
                        continue     
                    '''
                    if ( ( ( r1 in R1[pr][it]) or (R1[pr][it] in r1) ) \
                       and ( (r2 in R2[it]) or ( R2[it] in r2) ) ) \
                        or ( ( ( r1 in R2[it]) or (R2[it] in r1) ) \
                       and ( ( r2 in R1[pr][it]) or ( R1[pr][it] in r2 ) ) ):
                    '''
                    if (  ( r1 in R1[pr][it]) or (R1[pr][it] in r1)  \
                       
                        or ( r1 in R2[it]) or (R2[it] in r1)  ):
                        
                        
                        if r2 in R2[it] or R2[it] in r2 or r2 in R1[pr][it] or R1[pr][it] in r2 :
                            
                            if  ( ( r1.startswith(R1[pr][it]) or R1[pr][it].startswith(r1) ) \
                               and ( r2.startswith(R2[it]) or R2[it].startswith(r2) ) ) \
                                or ( ( r1.startswith(R2[it]) or R2[it].startswith(r1) ) \
                               and ( r2.startswith(R1[pr][it]) or R1[pr][it].startswith(r2)) ):
                            
                                
                                duplicates_local.append(it)
                    else:
                        break
                        
                        #pass
            print "Duplicates number = " ,
            
            tmp = set(duplicates_local)
            print len(tmp) ,
            count_dup += len(set(duplicates_local))
            print "Total duplicates  = " ,
            print count_dup
            
            duplicates.extend(tmp)
            duplicates_by_prefix[pr] = tmp
            del sorted_reads
            
            '''
                if it in duplicate_pairs and id in duplicate_pairs[it]:
                    continue
                
                if id not in duplicate_pairs:
                    duplicate_pairs[id]=[]
                    
                duplicate_pairs[id].append(it)  
            ''' 
                    
        
        print "total number of read pairs " ,
        print len(R2.keys()) ,
        print "total number of duplicates " ,
        print len(set(duplicates))
        
        with open(fl_duplicates, "w+") as outfile:             
            for pr in duplicates_by_prefix:
                for it in duplicates_by_prefix[pr]:
                    outfile.write(pr + "\t" + it + "\n")
        '''            
        with open(fl_R1, 'r') as data, \
             open(fl_R1_rmdup, 'w') as outfile:
            for lines in itertools.izip_longest(*[data]*4):  
                id=lines[0].rstrip().split(" ")[0]
                              
                pr=lines[1].rstrip()[0:9] 
                if id not in duplicates[pr]:
                    outfile.write(lines[0])
                    outfile.write(lines[1])
                    outfile.write(lines[2])
                    outfile.write(lines[3])
                    
        with open(fl_R2, 'r') as data, \
             open(fl_R2_rmdup, 'w') as outfile:
            for lines in itertools.izip_longest(*[data]*4):  
                id=lines[0].rstrip().split(" ")[0]
                pr=lines[1].rstrip()[0:9] 
                
                if id not in duplicates[pr]:
                    outfile.write(lines[0])
                    outfile.write(lines[1])
                    outfile.write(lines[2])
                    outfile.write(lines[3])        
        '''
                        
    #

#
#
def make_R1_compliment():
    folder="hic/raw_files_nongz/"
    
    filenames = [j for j in os.listdir(folder) if j.endswith(".fastq") and "_R1" in j]
    for fl in filenames:
        outfl=fl.split(".fastq")[0] + "_comp.fastq"
        with open(folder + fl, 'r') as data,\
             open (folder + outfl, 'w') as outdata:
            for lines in itertools.izip_longest(*[data]*4):
                
                read = returnCompliment(lines[1].rstrip())
                qual = lines[3].rstrip()
                outdata.write(lines[0])
                outdata.write(read + "\n")
                outdata.write(lines[2])
                outdata.write(qual + "\n")
                
                
#
#
def get_max_len_for_pair():
    barcodes_idx={}
    with open("barcodes_indexed", 'r') as bcs:
        for i, l in enumerate(bcs):
            l = l.rstrip().split("\t")
            barcodes_idx[int(l[0])] = l[1]

    for idx in barcodes_idx:  
        if idx % 2 != 0:
            continue        
                
        if idx/2 not in [35, 37, 50, 58, 60, 85, 88, 89]:
            continue
        
        fl_readsR1 = "hic/raw_files_nongz/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_reads_merged_R1_comp.fastq" 
        fl_readsR2 = "hic/raw_files_nongz/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_reads_merged_R2.fastq"         
        
        maxlen=0        
        with open(fl_readsR1, 'r') as data:
            for lines in itertools.izip_longest(*[data]*4):                
                read=lines[1].rstrip()
                if len(read) > maxlen:
                    maxlen=len(read)
        with open(fl_readsR2, 'r') as data:
            for lines in itertools.izip_longest(*[data]*4):                
                read=lines[1].rstrip()
                if len(read) > maxlen:
                    maxlen=len(read)        
       
        print 'idx = ' ,             
        print idx/2  ,  
        print " max len = " ,
        print maxlen
                

#
#
def make_same_length(maxlen):
    barcodes_idx = {}
    with open("barcodes_indexed", 'r') as bcs:
        for i, l in enumerate(bcs):
            l = l.rstrip().split("\t")
            barcodes_idx[int(l[0])] = l[1]

    for idx in barcodes_idx:  
        if idx % 2 != 0:
            continue        
                
        if idx/2 not in [35, 37, 50, 58, 60, 85, 88, 89]:
            continue
        
        fl_readsR1 = "hic/raw_files_nongz/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_reads_merged_R1_comp.fastq" 
        fl_readsR2 = "hic/raw_files_nongz/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_reads_merged_R2.fastq"         
        
        fl_readsR1_ajc = "hic/raw_files_nongz/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_reads_merged_len_adj_R1.fastq" 
        fl_readsR2_ajc = "hic/raw_files_nongz/" + str(idx/2) + "_" + barcodes_idx[idx][0:8] + "_reads_merged_len_adj_R2.fastq"   
        
               
        with open(fl_readsR1, 'r') as data, \
             open (fl_readsR1_ajc, 'w') as output:
            for lines in itertools.izip_longest(*[data]*4):                
                read = lines[1].rstrip()
                qual = lines[3].rstrip()                
                if len(read) < 15:
                    continue
                #print lines[0]
                #print read
                #print qual                
                for i in range(len(read), maxlen):
                    read = read + "N"
                    qual = qual + "6" 
                #print read
                #print qual
                output.write(lines[0])
                output.write(read + "\n")
                output.write(lines[2])
                output.write(qual + "\n")
                
                
        with open(fl_readsR2, 'r') as data, \
             open (fl_readsR2_ajc, 'w') as output:
            for lines in itertools.izip_longest(*[data]*4):                
                read = lines[1].rstrip()
                qual = lines[3].rstrip()                
                if len(read) < 15:
                    continue
                #print read
                #print qual
                for i in range(len(read), maxlen):
                    read = read + "N"
                    qual = qual + "6" 
                #print read
                #print qual
                output.write(lines[0])
                output.write(read + "\n")
                output.write(lines[2])
                output.write(qual + "\n")        
       

                    
#
#
'''
create_automaton()
(B_fw, B_rv) = load_automaton_from_pickle()

for key in B_fw:
    idx = B_fw.get(key)[0]
    print  idx,
    print key ,    
    print B_rv.get(returnReverseCompliment(key))[0] ,
    print returnReverseCompliment(key)
'''

#folder="Ulianov-1-13_R1_001.fastq.gz"
#sorting_with_automaton(folder, B_fw, B_rv)

#folder="Ulianov-1-13_R2_001.fastq.gz"
#sorting_with_automaton(folder, B_fw, B_rv)


# zcat Ulianov-1-13_R2_001.fastq.gz | python ahocorasick_search.py 
#download_from_s3()


#parse_reads()
#concatenate_reads()
#check_duplicates()

#make_R1_compliment()
#get_max_len_for_pair()

#make_same_length(117)
'''
Created on May 12, 2017
@author: ky6, pika1
'''
import os, datetime, re, csv
import pandas as pd
from openpyxl import Workbook
from openpyxl import load_workbook
from collections import Counter
import numpy as np
import xlrd

####
# change input file directory here, and change the format and order of drugs
####
# ranks = {"mod*nsd*opd":4,"nsd*opd":4,"mod*opd":3,"opd":3,"mod*nsd":2,"nsd":2,"mod":1}
# ranks_bigram = {"opd-opd":2,"nsd-nsd-opd-nsd-opd-nsd-mod":1,"mod-nsd-opd":4,"nsd-mod-opd":4,"opd-nsd-mod":4,"nsd-opd":4,"opd-nsd":4,"mod-opd":3,"opd":3,"mod-nsd":2,"nsd-mod":2,"nsd":2,"mod":1,"opd-mod":3,"opd-mod-opd":3,"nsd-mod-nsd-mod-opd-nsd":1}
# # ranks_bigram = {"opd-opd":2}

base_dir = "~/gra"
file_name = 'N_O_M_seq.csv'
first_occurence = {}
last_occurence = {}
first_useage = {}
last_useage = {}

drugs = ["mod","nsd","opd","mod*nsd","mod*opd","nsd*opd","mod*nsd*opd"]
bigram_patterns = []
bigram = {}
bigram_row_frequency = {}
bigram_frequency = {}
bigram_row_percent = {}
bigram_percent = {}
bigram_min = {}
bigram_max = {}
bigram_weighted_avg = {}
bigram_avg = {}
ranks ={}
bigram = {}
bigram_pair_freq_array ={}
row_pattern_count = {}
bigram_pair_weighted_array = {}
mean_array = {}
square_array = {}
pattern_count = {}
raw_row_count = {}
flag_pattern = {}

def Read_Excel():
    '''
    load data from excel using pandas
    '''
    with open(file_name,'rb') as csvfile:
        wsheet = csv.reader(csvfile,delimiter=',')
        mydata = []
        total_count = np.float(0)
        for index, row in enumerate(wsheet):
            drug_seq = row[0]
            freq = int(row[1])
            percent = row[2]
            # freq = wsheet['Frequency'][row]
            total_count += int(freq)
            # percent = wsheet['Percent'][row]
            mydata.append({'drug_seq':drug_seq,'freq':freq,'percent':percent})
    return mydata, total_count

def Insert_Opd(w):
    idx = w.find('op')
    if idx==-1:
        return w
    if idx+2<=len(w)-1 and w[idx+2]=='d':
        return w
    else:
        return w[:idx+2]+'d'+w[idx+2:]

def process_values(pattern,seq,d,s,type):
    rowcount = 0
    pattern_count[pattern] = 0
    if((pattern in bigram_patterns)):
        rowcount = rowcount + 1
        count = bigram[pattern] + 1
        row_pattern_count[pattern] = row_pattern_count[pattern] + 1
        bigram[pattern] = count

        if(bigram_row_frequency[pattern] == 0):
            bigram_frequency[pattern] = d['freq'] + bigram_frequency[pattern]

        if(bigram_row_percent[pattern] == 0):
            bigram_percent[pattern] = float(d['percent']) + float(bigram_percent[pattern])

        bigram_row_frequency[pattern] = d['freq']
        bigram_row_percent[pattern] = d['percent']
         
        if(bigram_min[pattern] == 0):
            bigram_min[pattern] = row_pattern_count[pattern]

        if(bigram_max[pattern] < row_pattern_count[pattern]):
            bigram_max[pattern] = row_pattern_count[pattern]

        #first occurence index
        if pattern in first_occurence and (s == 0):
            first_occurence[pattern] = 1 + first_occurence[pattern]
            first_useage[pattern] = first_useage[pattern] + d['freq']

        #last occurence index
        if(type == "monogram"):
            param = (len(seq)-1)
        elif(type == "bigram"):
            param = (len(seq)-2)
        elif(type == "trigram"):
            param = (len(seq)-3)

        if(s == param):
            if pattern in last_occurence:
                last_occurence[pattern] = 1 + last_occurence[pattern]
                last_useage[pattern] = last_useage[pattern] + d['freq']

        endval = 0
        if(s == param):
            if(type == "monogram"):
                endval = len(seq)
            elif(type == "bigram"):
                endval = len(seq) - 1
            elif(type == "trigram"):
                endval = len(seq) -2
            
            for r in range(0,endval):
                if(type == "monogram"):
                    pattern  = (seq[r])
                elif(type == "bigram"):
                    pattern = (seq[r]+'-'+seq[r+1])
                elif(type == "trigram"):
                    pattern = (seq[r]+'-'+seq[r+1]+'-'+seq[r+2])

                # print pattern,pattern_count[pattern]
                if pattern_count[pattern] == (row_pattern_count[pattern] - 1):
                    # print pattern, row_pattern_count[pattern], bigram_row_frequency[pattern]
                    bigram_pair_weighted_array[pattern].append((float(bigram_row_frequency[pattern] * row_pattern_count[pattern])))
                    raw_row_count[pattern].append(row_pattern_count[pattern])
                else:
                    pattern_count[pattern] = pattern_count[pattern] + 1

    if(flag_pattern[pattern] == 0):
        bigram_pair_freq_array[pattern].append(d['freq'])
        flag_pattern[pattern] = 1

def calculate_std_dev(type):
    pattern_array = []
    if(type == "bigram" or type == "trigram"):
        pattern_array = bigram_patterns
    elif(type == "monogram"):
        pattern_array = drugs

    for index, d in enumerate(pattern_array):
        if(bigram[pattern_array[index]] == 0):
            bigram_pair_weighted_array[pattern_array[index]] = 0
        else:
            for j in range(0,len(bigram_pair_weighted_array[pattern_array[index]])):
                mean_array[pattern_array[index]] = mean_array[pattern_array[index]] + bigram_pair_weighted_array[pattern_array[index]][j]
            mean_array[pattern_array[index]] = round((mean_array[pattern_array[index]] / bigram[pattern_array[index]]),2)

    square = 0
    for index, d in enumerate(pattern_array):
        square_array[pattern_array[index]] = 0
        if(bigram[pattern_array[index]] == 0):
            bigram_avg[pattern_array[index]] = 0
        else:
            for j in range(0,len(raw_row_count[pattern_array[index]])):
                # print pattern_array[index],raw_row_count[pattern_array[index]],mean_array[pattern_array[index]],bigram_pair_freq_array[pattern_array[index]]
                square = (((raw_row_count[pattern_array[index]][j] - mean_array[pattern_array[index]]) * (raw_row_count[pattern_array[index]][j] - mean_array[pattern_array[index]])) * bigram_pair_freq_array[pattern_array[index]][j])
                square_array[pattern_array[index]] = square_array[pattern_array[index]] + square
                # print "square is ", square
                # print "squre sum", square_array[pattern_array[index]]
        
            # print bigram_pair_weighted_array[pattern_array[index]],mean_array[pattern_array[index]]
            
            # print "squre sum", square_array[pattern_array[index]],"count",bigram[pattern_array[index]]
            bigram_avg[pattern_array[index]] = (square_array[pattern_array[index]] / bigram[pattern_array[index]])
            # print "divide is ",pattern_array[index], bigram_avg[pattern_array[index]]
            bigram_avg[pattern_array[index]] = (bigram_avg[pattern_array[index]])**(1/2.0)
            # print "std dev",pattern_array[index], bigram_avg[pattern_array[index]]


def clearRowCount():
    for i in range(0, len(drugs)):
        for j in range(0,len(drugs)):
            pattern = (drugs[i]+'-'+drugs[j])
            row_pattern_count[pattern] = 0
            bigram_row_frequency[pattern] = 0
            bigram_row_percent[pattern] = 0
            flag_pattern[pattern] = 0

def clearMonoRowCount():
    for i in range(0, len(drugs)):
        row_pattern_count[(drugs[i])] = 0
        bigram_row_frequency[(drugs[i])] = 0
        bigram_row_percent[(drugs[i])] = 0
        pattern_count[(drugs[i])] = 0
        flag_pattern[(drugs[i])] = 0

def clearTrigramRowCount():
    for i in range(0, len(drugs)):
        for j in range(0,len(drugs)):
            for k in range(0,len(drugs)):
                row_pattern_count[(drugs[i]+'-'+drugs[j]+'-'+drugs[k])] = 0
                bigram_row_frequency[(drugs[i]+'-'+drugs[j]+'-'+drugs[k])] = 0
                bigram_row_percent[(drugs[i]+'-'+drugs[j]+'-'+drugs[k])] = 0
                flag_pattern[(drugs[i]+'-'+drugs[j]+'-'+drugs[k])] = 0

def Process(data,total_count):
    Pattern_generator()
    type = "monogram"
    result = {}
    rowcount = 0
    for index, d in enumerate(data):
        seq = [Insert_Opd(w) for w in d['drug_seq'].split('-')]
        clearMonoRowCount()
        for s in range(0,len(seq)):
            if(s != len(seq)):
                pattern  = (seq[s])
                process_values(pattern,seq,d,s,type)
        
    calculate_std_dev("monogram")
                
    aggregate_value = open('monogram.csv','wb')
    print >> aggregate_value, "Index," , "Pattern," ,"Count,", "Frequency,", "Percent,", "Min,","Max,","Weighted average,","Std. Dev.,", "First occurrence count,","Last occurence count,","First occurence,","Last occurence"
    rowcount = 0

    for index, data in enumerate(drugs):
        if(bigram[drugs[index]] != 0):
            rowcount = rowcount +1
            print >> aggregate_value, rowcount,',', drugs[index],',',bigram[drugs[index]],',', bigram_frequency[drugs[index]],',', bigram_percent[drugs[index]],',',bigram_min[drugs[index]],',',bigram_max[drugs[index]],',',round(mean_array[drugs[index]],2),',',round(bigram_avg[drugs[index]],2),',',first_occurence[drugs[index]],',',last_occurence[drugs[index]],',',first_useage[drugs[index]],',',last_useage[drugs[index]]
    aggregate_value.close()

def Process_bigram(data,total_count):
    generate_bigram()
    type = "bigram"
    result = {}
    max_mapping = {}
    for index, d in enumerate(data):
        seq = [Insert_Opd(w) for w in d['drug_seq'].split('-')]
        clearRowCount()
        for s in range(0,len(seq)-1):
            if(s != len(seq)):
                pattern  = (seq[s]+'-'+(seq[s+1]))
                process_values(pattern,seq,d,s,type)

    calculate_std_dev(type)

    aggregate_value = open('bigram.csv','wb')
    print >> aggregate_value, "Index," , "Pattern," ,"Count,", "Frequency,", "Percent,", "Min,","Max,","Weighted average,","Std. Dev.,", "First occurrence count,","Last occurence count,","First occurence,","Last occurence"
    rowcount = 0

    for index, data in enumerate(bigram_patterns):
        if(bigram[bigram_patterns[index]] != 0 and (bigram_patterns[index].count("-") == 1)):
            rowcount = rowcount +1
            print >> aggregate_value, rowcount,',', bigram_patterns[index],',',bigram[bigram_patterns[index]],',', bigram_frequency[bigram_patterns[index]],',', bigram_percent[bigram_patterns[index]],',',bigram_min[bigram_patterns[index]],',',bigram_max[bigram_patterns[index]],',',round(mean_array[bigram_patterns[index]],2),',',round(bigram_avg[bigram_patterns[index]],2),',',first_occurence[bigram_patterns[index]],',',last_occurence[bigram_patterns[index]],',',first_useage[bigram_patterns[index]],',',last_useage[bigram_patterns[index]]
    aggregate_value.close()
        
def Process_Trigram(data,total_count):
    generate_trigram()
    type = "trigram"
    generate_trigram()
    result = {}
    for index, d in enumerate(data):
        seq = [Insert_Opd(w) for w in d['drug_seq'].split('-')]
        clearTrigramRowCount()
        for s in range(0,len(seq)-1):
            if(s != len(seq)-2):
                pattern  = (seq[s]+'-'+(seq[s+1])+'-'+seq[s+2])
                process_values(pattern,seq,d,s,type)

    calculate_std_dev(type)

    aggregate_value = open('trigram.csv','wb')
    print >> aggregate_value, "Index," , "Pattern," ,"Count,", "Frequency,", "Percent,", "Min,","Max,","Weighted average,","Std. Dev.,", "First occurrence count,","Last occurence count,","First occurence,","Last occurence"
    rowcount = 0

    for index, data in enumerate(bigram_patterns):
        if(bigram[bigram_patterns[index]] != 0 and (bigram_patterns[index].count("-") == 2)):
            rowcount = rowcount +1
            print >> aggregate_value, rowcount,',', bigram_patterns[index],',',bigram[bigram_patterns[index]],',', bigram_frequency[bigram_patterns[index]],',', bigram_percent[bigram_patterns[index]],',',bigram_min[bigram_patterns[index]],',',bigram_max[bigram_patterns[index]],',',round(mean_array[bigram_patterns[index]],2),',',round(bigram_avg[bigram_patterns[index]],2),',',first_occurence[bigram_patterns[index]],',',last_occurence[bigram_patterns[index]],',',first_useage[bigram_patterns[index]],',',last_useage[bigram_patterns[index]]
    aggregate_value.close()
    
def generate_trigram():
    for i in range(0, len(drugs)):
        for j in range(0,len(drugs)):
            for k in range(0,len(drugs)):
                pattern = (drugs[i]+'-'+drugs[j]+'-'+drugs[k])
                first_occurence[pattern] = 0
                last_occurence[pattern] = 0
                first_useage[pattern] = 0
                last_useage[pattern] = 0
                bigram[pattern] = 0
                bigram_frequency[pattern] = 0
                ranks[pattern] = 2
                bigram_patterns.append(pattern)
                bigram_row_frequency[pattern] = 0
                bigram_row_percent[pattern] = 0.0
                bigram_percent[pattern] = 0.0
                bigram_min[pattern] = 0
                bigram_max[pattern] = 0
                bigram_weighted_avg[pattern] = 0
                bigram_avg[pattern] = 0
                bigram_pair_freq_array[pattern] = []
                bigram_pair_weighted_array[pattern] = []
                row_pattern_count[pattern] = 0
                mean_array[pattern] = 0
                square_array[pattern] = 0
                raw_row_count[pattern] = []
                flag_pattern[pattern] = 0

def generate_bigram():
    for i in range(0, len(drugs)):
        for j in range(0,len(drugs)):
            pattern = (drugs[i]+'-'+drugs[j])
            first_occurence[pattern] = 0
            last_occurence[pattern] = 0
            first_useage[pattern] = 0
            last_useage[pattern] = 0
            bigram[pattern] = 0
            bigram_frequency[pattern] = 0
            ranks[pattern] = 2
            bigram_patterns.append(pattern)
            bigram_row_frequency[pattern] = 0
            bigram_row_percent[pattern] = 0.0
            bigram_percent[pattern] = 0.0
            bigram_min[pattern] = 0
            bigram_max[pattern] = 0
            bigram_weighted_avg[pattern] = 0
            bigram_avg[pattern] = 0
            bigram_pair_freq_array[pattern] = []
            bigram_pair_weighted_array[pattern] = []
            row_pattern_count[pattern] = 0
            mean_array[pattern] = 0
            square_array[pattern] = 0
            raw_row_count[pattern] = []
            flag_pattern[pattern] = 0


def Pattern_generator():
    for i in range(0, len(drugs)):
        pattern = drugs[i]
        first_occurence[pattern] = 0
        last_occurence[pattern] = 0
        first_useage[pattern] = 0
        last_useage[pattern] = 0
        bigram[pattern] = 0
        bigram_frequency[pattern] = 0
        ranks[pattern] = 2
        bigram_patterns.append(pattern)
        bigram_row_frequency[pattern] = 0
        bigram_row_percent[pattern] = 0.0
        bigram_percent[pattern] = 0.0
        bigram_min[pattern] = 0
        bigram_max[pattern] = 0
        bigram_weighted_avg[pattern] = 0
        bigram_avg[pattern] = 0
        bigram_pair_freq_array[pattern] = []
        bigram_pair_weighted_array[pattern] = []
        row_pattern_count[pattern] = 0
        mean_array[pattern] = 0
        square_array[pattern] = 0
        raw_row_count[pattern] = []
        flag_pattern[pattern] = 0

def Main():
    data,total_count = Read_Excel()
    print("data is:",data)
    print("total is:",total_count)
    Process(data,total_count)
    Process_bigram(data,total_count)
    Process_Trigram(data,total_count)
    

if __name__ == '__main__':
    Main()


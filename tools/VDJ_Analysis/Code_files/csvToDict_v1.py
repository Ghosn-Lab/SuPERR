# -*- coding: utf-8 -*-
"""
Spyder Editor

This script generates a dictionary of mutations associated with a cell
"""
import json
fname = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\VDJ\\superrseq_processed_detailed_mutation_data.csv"
fileCont = []
with open(fname, 'r') as seqH:
    fileCont = seqH.readlines()

contigDict = {}
for entry in fileCont[1:]:
    contDet = entry.split(",")
    contigId = contDet[0].split("_")[0]
    if contigId not in contigDict:
        contigDict[contigId] = []
        contigDict[contigId].append({'region': contDet[2].split("-")[0], 
                                     'position': contDet[3],
                                     'germline': contDet[4],
                                     'mutation': contDet[5],
                                     'isHS': contDet[6],
                                     'hsType': contDet[7],
                                     'isCoding': contDet[8]})
    else:
        contigDict[contigId].append({'region': contDet[2].split("-")[0], 
                                     'position': contDet[3],
                                     'germline': contDet[4],
                                     'mutation': contDet[5],
                                     'isHS': contDet[6],
                                     'hsType': contDet[7],
                                     'isCoding': contDet[8]})

wfName = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\VDJ\\dictionary.txt"
with open(wfName, 'w') as wFile:
    wFile.write(json.dumps(contigDict))
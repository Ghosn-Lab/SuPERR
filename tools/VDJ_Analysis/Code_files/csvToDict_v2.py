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

regionWiseInfo = {}
hsTypeInfo = {}
for cell, muts in contigDict.items():
    for mut in muts:
        if mut["region"] not in regionWiseInfo and mut["region"] not in hsTypeInfo:
            hsTypeInfo[mut["region"]] = {}
            if mut["hsType"] != "NA":
                if mut["hsType"] not in hsTypeInfo[mut["region"]]:
                    hsTypeInfo[mut["region"]][mut["hsType"]] = 1
                else:
                    hsTypeInfo[mut["region"]][mut["hsType"]] = hsTypeInfo[mut["region"]][mut["hsType"]] + 1
            regionWiseInfo[mut["region"]] = {}
            hsCount = (0, 1) [mut["isHS"] == "TRUE"]
            regionWiseInfo[mut["region"]][cell] = [1, hsCount]
        else:
            if mut["hsType"] != "NA":
                if mut["hsType"] not in hsTypeInfo[mut["region"]]:
                    hsTypeInfo[mut["region"]][mut["hsType"]] = 1
                else:
                    hsTypeInfo[mut["region"]][mut["hsType"]] = hsTypeInfo[mut["region"]][mut["hsType"]] + 1
            if cell not in regionWiseInfo[mut["region"]]:
                hsCount = (0, 1) [mut["isHS"] == "TRUE"]
                regionWiseInfo[mut["region"]][cell] = [1, hsCount]
            else:
                mutCount = regionWiseInfo[mut["region"]][cell][0] + 1
                hsCount = regionWiseInfo[mut["region"]][cell][1]
                hsCount = (hsCount, hsCount+1) [mut["isHS"] == "TRUE"]
                regionWiseInfo[mut["region"]][cell] = [mutCount, hsCount]


# for cell, muts in contigDict.items():
#     for mut in muts:
#         if mut["region"] not in hsTypeInfo:
            
#         else:
            

wfName0 = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\VDJ\\contigs_data.json"
with open(wfName0, 'w') as wFile:
    wFile.write(json.dumps(contigDict))
    
wfName1 = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\VDJ\\region_data.json"
with open(wfName1, 'w') as wFile:
    wFile.write(json.dumps(regionWiseInfo))
    
wfName2 = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\VDJ\\hs_data.json"
with open(wfName2, 'w') as wFile:
    wFile.write(json.dumps(hsTypeInfo))
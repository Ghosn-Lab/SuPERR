# -*- coding: utf-8 -*-
"""
Spyder Editor

This script generates a dictionary of mutations associated with a cell
"""
import json, csv
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
hsTypeInfo = {'name': 'regions', 'children': []}
for cell, muts in contigDict.items():
    for mut in muts:
        if mut["region"] not in regionWiseInfo:
            regionWiseInfo[mut["region"]] = []
            hsCount = (0, 1) [mut["isHS"] == "TRUE"]
            regionWiseInfo[mut["region"]].append({'name': cell, 'total_muts': 1, 'hs_muts': hsCount})
        else:
            isCell = list(cellN for cellN in regionWiseInfo[mut["region"]] if cellN['name'] == cell)
            if isCell == []:
                hsCount = (0, 1) [mut["isHS"] == "TRUE"]
                regionWiseInfo[mut["region"]].append({'name': cell, 'total_muts': 1, 'hs_muts': hsCount})
            else:
                isCell[0]['total_muts'] += 1
                isCell[0]['hs_muts'] = (isCell[0]['hs_muts'], isCell[0]['hs_muts'] + 1) [mut["isHS"] == "TRUE"]


for cell, muts in contigDict.items():
    for mut in muts:
        isRegion = list(region for region in hsTypeInfo['children'] if region['name'] == mut["region"])
        if isRegion == []:
            hsTypeInfo['children'].append({'name': mut["region"], 'children': []})
            if mut["hsType"] != "NA":
                region = list(region for region in hsTypeInfo['children'] if region['name'] == mut["region"])
                isHS = list(hs for hs in region[0]['children'] if hs['name'] == mut["hsType"])
                if isHS == []:
                    region[0]['children'].append({'name': mut["hsType"], 'value': 1})
                else:
                    isHS[0]['value'] = isHS[0]['value'] + 1
        else:
            if mut["hsType"] != "NA":
                region = list(region for region in hsTypeInfo['children'] if region['name'] == mut["region"])
                isHS = list(hs for hs in region[0]['children'] if hs['name'] == mut["hsType"])
                if isHS == []:
                    region[0]['children'].append({'name': mut["hsType"], 'value': 1})
                else:
                    isHS[0]['value'] = isHS[0]['value'] + 1
            

# wfName0 = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\VDJ\\contigs_data.json"
# with open(wfName0, 'w') as wFile:
#     wFile.write(json.dumps(contigDict))
    
wfName1 = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\VDJ\\region_data.json"
with open(wfName1, 'w') as wFile:
    wFile.write(json.dumps(regionWiseInfo))

def writeToCSV(rFName, listDict):
    headers = listDict[0].keys()
    with open(rFName, 'w', newline = '') as outFile:
        csvDictWriter = csv.DictWriter(outFile, headers)
        csvDictWriter.writeheader()
        csvDictWriter.writerows(listDict)

wfNameCSV = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\VDJ\\"
for region in regionWiseInfo:
    regionCellList = regionWiseInfo[region]
    outFName = wfNameCSV + region + ".csv"
    writeToCSV(outFName, regionCellList)
    
wfName2 = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\VDJ\\hs_data.json"
with open(wfName2, 'w') as wFile:
    wFile.write(json.dumps(hsTypeInfo))
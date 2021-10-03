# -*- coding: utf-8 -*-
"""
Spyder Editor

This script generates a dictionary of mutations associated with a cell
"""
import json, csv, os
# processes the detailed mutation data for each barcode
def mutationDataProcessing(inFileNames, outDir):
    fname = inFileNames[0]
    fileCont = []
    with open(fname, 'r') as seqH:
        fileCont = seqH.readlines()
    # get the mutation related to each contig
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
    # get the total mutation data from the base imgt output
    fname1 = inFileNames[1]
    totalFileCont = []
    with open(fname1, 'r') as seqH:
        totalFileCont = seqH.readlines()
    # generate cumulative results for mutations for the contigs
    allContigDict = {}
    for entry in totalFileCont[1:]:
        contDet = entry.split(",")
        contigId = contDet[0].split("_")[0]
        if contigId not in allContigDict:
            allContigDict[contigId] = {}
            allContigDict[contigId] = {'FR1': int(contDet[377]), 
                                       'FR2': int(contDet[378]),
                                       'FR3': int(contDet[379]),
                                       'CDR1': int(contDet[380]),
                                       'CDR2': int(contDet[381]),
                                       'CDR3': int(contDet[382])}
        else:
            allContigDict[contigId]['FR1'] += int(contDet[377]) 
            allContigDict[contigId]['FR2'] += int(contDet[378])
            allContigDict[contigId]['FR3'] += int(contDet[379])
            allContigDict[contigId]['CDR1'] += int(contDet[380])
            allContigDict[contigId]['CDR2'] += int(contDet[381])
            allContigDict[contigId]['CDR3'] += int(contDet[382])
    
    # create region wise mutation info (CDR or FR)
    regionWiseInfo = {}
    hsTypeInfo = {'name': 'regions', 'children': []}
    for cell, muts in contigDict.items():
        for mut in muts:
            if mut["region"] not in regionWiseInfo:
                regionWiseInfo[mut["region"]] = []
                hsCount = (0, 1) [mut["isHS"] == "TRUE"]
                regionWiseInfo[mut["region"]].append({'name': cell, 'total_muts': 1, 'hs_muts': hsCount, 'nb_of_hotspots': allContigDict[cell][mut["region"]]})
            else:
                isCell = list(cellN for cellN in regionWiseInfo[mut["region"]] if cellN['name'] == cell)
                if isCell == []:
                    hsCount = (0, 1) [mut["isHS"] == "TRUE"]
                    regionWiseInfo[mut["region"]].append({'name': cell, 'total_muts': 1, 'hs_muts': hsCount, 'nb_of_hotspots': allContigDict[cell][mut["region"]]})
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
                
    
    # wfName0 = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\VDJ\\Results_Prelims\\contigs_data.json"
    # with open(wfName0, 'w') as wFile:
    #     wFile.write(json.dumps(allContigDict))
    # write to region wise mutation json output    
    wfName1 = os.path.join(outDir, "region_data.json")
    with open(wfName1, 'w') as wFile:
        wFile.write(json.dumps(regionWiseInfo))
    
    def writeToCSV(rFName, listDict):
        headers = listDict[0].keys()
        with open(rFName, 'w', newline = '') as outFile:
            csvDictWriter = csv.DictWriter(outFile, headers)
            csvDictWriter.writeheader()
            csvDictWriter.writerows(listDict)
    # outputs region wise mutation data - cumulative
    mutFilePaths = []
    wfNameCSV = os.path.join(outDir, "temp")
    for region in regionWiseInfo:
        regionCellList = regionWiseInfo[region]
        outFName = region + ".csv"
        outFilePath = os.path.join(wfNameCSV, outFName)
        mutFilePaths.append(outFilePath)
        writeToCSV(outFilePath, regionCellList)
    # write to hotspot details region wise    
    wfName2 = os.path.join(outDir, "Mutation_Data", "data", "hs_data.json")
    with open(wfName2, 'w') as wFile:
        wFile.write(json.dumps(hsTypeInfo))
    
    return mutFilePaths
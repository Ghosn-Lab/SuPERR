# -*- coding: utf-8 -*-
"""
@author: scsac
"""

import csv, json, os
# generates a circos plotting json file from imgt output
def circosFileGenerator(inFileName, outDir):
    tempFileDir = os.path.join(outDir, "temp")
    tempFileName = os.path.join(tempFileDir, "lineage_processed_temp.csv")
    # create temp file to remove repeat contig barcodes
    tempFileContent = []
    with open(inFileName, 'r') as tempR:
        tempFileContent = tempR.readlines()
    # headers required are : barcode, isotype and lineage-id
    headers = ["indexed_barcode", "isotype", "lineage_id"]
    headerInd = []
    headerList = tempFileContent[0].rstrip().split(',')
    for header in headers:
        headerInd.append(headerList.index(header))
    
    indBarcodeData = []
    for barcode in tempFileContent[1:]:
        barCont = barcode.rstrip().split(',')
        isBarcode = list(bC for bC in indBarcodeData if bC[headers[0]] == barCont[headerInd[0]])
        if isBarcode == []:
            indBarcodeData.append({headers[0]: barCont[headerInd[0]], headers[1]: barCont[headerInd[1]], headers[2]: barCont[headerInd[2]]})
    
    with open(tempFileName, 'w', newline = '') as tempOut:
        csvDictWriter = csv.DictWriter(tempOut, headers)
        csvDictWriter.writeheader()
        csvDictWriter.writerows(indBarcodeData)
    
    filename = tempFileName
    # use the temp file to generate the circos matrix for the data
    isotypeGroups = set()
    jsonContentUnSorted = {}
    
    with open(filename, 'r') as csvH:
        csvReader = csv.reader(csvH)
        for lineageData in csvReader:
            if lineageData[1] not in jsonContentUnSorted:
                jsonContentUnSorted[lineageData[1]] = []
                jsonContentUnSorted[lineageData[1]].append(lineageData[2])
            else:
                jsonContentUnSorted[lineageData[1]].append(lineageData[2])
            isotypeGroups.add(lineageData[1])
    # gather isotype information and groups
    isotypeGroups.remove('isotype')
    isotypeGroups = sorted(list(isotypeGroups))
    
    jsonContentSorted = []
    for isotype in isotypeGroups:
        jsonContentSorted.append({'name': isotype, 'matrix': [0]*len(isotypeGroups)})
    # create matrix with shared lineages across isotypes
    sharedLineageIDs = {}
    i = 0
    while i < len(isotypeGroups):
        for isotype in isotypeGroups[i+1:]:
            sharedLI = set(jsonContentUnSorted[isotypeGroups[i]]).intersection(set(jsonContentUnSorted[isotype]))
            if len(sharedLI) != 0:
                if isotypeGroups[i] not in sharedLineageIDs:
                    sharedLineageIDs[isotypeGroups[i]] = list(sharedLI)
                else:
                    sharedLineageIDs[isotypeGroups[i]] += list(sharedLI)
                iso1count = 0
                iso2count = 0
                for lineage in sharedLI:
                    iso1count += jsonContentUnSorted[isotypeGroups[i]].count(lineage)
                    iso2count += jsonContentUnSorted[isotype].count(lineage)
                isoItem1 = list(isoT['matrix'] for isoT in jsonContentSorted if isoT['name'] == isotypeGroups[i])
                isoItem2 = list(isoT['matrix'] for isoT in jsonContentSorted if isoT['name'] == isotype)
                isoItem1[0][isotypeGroups.index(isotype)] = iso1count
                isoItem2[0][isotypeGroups.index(isotypeGroups[i])] = iso2count
        i += 1
    # refine matrix with shared lineage ids
    for isotype in isotypeGroups:
        lineageCount = 0
        if isotype in sharedLineageIDs:
            totalNonSharedLIds = set(jsonContentUnSorted[isotype]) - set(sharedLineageIDs[isotype])
            lineageCount = len(totalNonSharedLIds)
        else:
            lineageCount = len(set(jsonContentUnSorted[isotype]))
        isoItem = list(isoT['matrix'] for isoT in jsonContentSorted if isoT['name'] == isotype)
        isoItem[0][isotypeGroups.index(isotype)] = lineageCount
    # write json file with matrix data to output
    with open(os.path.join(outDir, "Circos_Plot", "data", "circos.json"), 'w') as wFile:
        wFile.write(json.dumps(jsonContentSorted))
# -*- coding: utf-8 -*-
"""
@author: scsac
"""

import csv, os
from pathlib import Path

# generates the comparative files for VDJ data for samples included
def vdj_CompareFiles(inputFileNames, outDir):
    v_data = []
    d_data = []
    j_data = []
    header = ['name']
    for fileName in inputFileNames:
        header.append(os.path.basename(str(Path(fileName).parent)))
        fileContent = []
        with open(fileName, 'r') as contig_file:
            fileContent = contig_file.readlines()
        for barcode in fileContent:
            barCont = barcode.split(",")
            if barCont[11] == 'True':
                # check for v gene in v_data
                isV = list(vName for vName in v_data if vName['name'] == barCont[6])
                if isV == []:
                    v_data.append({'name': barCont[6], header[-1]: 1})
                else:
                    if header[-1] in isV[0]:
                        isV[0][header[-1]] = isV[0][header[-1]] + 1
                    else:
                        isV[0][header[-1]] = 1
                # check for d gene in d_data
                isD = list(dName for dName in d_data if dName['name'] == barCont[7])
                if isD == []:
                    d_data.append({'name': barCont[7], header[-1]: 1})
                else:
                    if header[-1] in isD[0]:
                        isD[0][header[-1]] = isD[0][header[-1]] + 1
                    else:
                        isD[0][header[-1]] = 1
                # check for j gene in j_data     
                isJ = list(jName for jName in j_data if jName['name'] == barCont[8])
                if isJ == []:
                    j_data.append({'name': barCont[8], header[-1]: 1})
                else:
                    if header[-1] in isJ[0]:
                        isJ[0][header[-1]] = isJ[0][header[-1]] + 1
                    else:
                        isJ[0][header[-1]] = 1
    
    # remove none from the data lists
    for i in range(len(v_data)):
        if v_data[i]['name'] == 'None':
            del v_data[i]
            break
    for i in range(len(d_data)):
        if d_data[i]['name'] == 'None':
            del d_data[i]
            break
    for i in range(len(j_data)):
        if j_data[i]['name'] == 'None':
            del j_data[i]
            break
    
    # write zeroes to files without any data for particular genes
    for file in header[1:]:
        for v_gene in v_data:
            if file not in v_gene:
                v_gene[file] = 0
        for d_gene in d_data:
            if file not in d_gene:
                d_gene[file] = 0
        for j_gene in j_data:
            if file not in j_gene:
                j_gene[file] = 0
    
    # write files to output directory
    wfName1 = os.path.join(outDir, "VDJ_Comparative", "data", "v_CompData.csv")
    with open(wfName1, 'w', newline="") as writeFile:
        writer = csv.DictWriter(writeFile, fieldnames=header)
        writer.writeheader()
        for v_gene in v_data:
            writer.writerow(v_gene)
    
    wfName2 = os.path.join(outDir, "VDJ_Comparative", "data", "d_CompData.csv")
    with open(wfName2, 'w', newline="") as writeFile:
        writer = csv.DictWriter(writeFile, fieldnames=header)
        writer.writeheader()
        for d_gene in d_data:
            writer.writerow(d_gene)
    
    wfName3 = os.path.join(outDir, "VDJ_Comparative", "data", "j_CompData.csv")
    with open(wfName3, 'w', newline="") as writeFile:
        writer = csv.DictWriter(writeFile, fieldnames=header)
        writer.writeheader()
        for j_gene in j_data:
            writer.writerow(j_gene)
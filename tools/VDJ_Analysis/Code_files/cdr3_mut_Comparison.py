# -*- coding: utf-8 -*-
"""
@author: scsac
"""

import csv, os
from pathlib import Path

# generates the comparative files for VDJ data for samples included
def mut_CompareFiles(inputFileNames, outDir):
    mut_data = []
    header = ['name']
    for fileName in inputFileNames:
        header.append(os.path.basename(str(Path(fileName).parent)))
        fileContent = []
        with open(fileName, 'r') as contig_file:
            fileContent = contig_file.readlines()
        for barcode in fileContent:
            barCont = barcode.split(",")
            if barCont[6] == 'TRUE':
                # check for cdr3 hotspot in mutation
                if barCont[2].split('-')[0] == 'CDR3':
                    isHS = list(hsName for hsName in mut_data if hsName['name'] == barCont[7])
                    if isHS == []:
                        mut_data.append({'name': barCont[7], header[-1]: 1})
                    else:
                        if header[-1] in isHS[0]:
                            isHS[0][header[-1]] = isHS[0][header[-1]] + 1
                        else:
                            isHS[0][header[-1]] = 1
    
    # remove none from the data lists
    for i in range(len(mut_data)):
        if mut_data[i]['name'] == 'None':
            del mut_data[i]
            break 
    
    # write zeroes to files without any data for particular genes
    for file in header[1:]:
        for mutHS in mut_data:
            if file not in mutHS:
                mutHS[file] = 0
    
    # write files to output directory
    wfName = os.path.join(outDir, "Mut_Comparative", "data", "mut_CompData.csv")
    with open(wfName, 'w', newline="") as writeFile:
        writer = csv.DictWriter(writeFile, fieldnames=header)
        writer.writeheader()
        for mutHS in mut_data:
            writer.writerow(mutHS)
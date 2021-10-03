# -*- coding: utf-8 -*-
"""
@author: scsac
"""

import csv, os, sys

# generates the comparative files for VDJ data for samples included
def mut_CompareBarcodes(fileName, outDir):
    mut_data = []
    header = ['name']
    fileContent = []
    with open(fileName, 'r') as contig_file:
        fileContent = contig_file.readlines()
    for barcode in fileContent:
        barCont = barcode.split(",")
        if barCont[6] == 'TRUE':
            # check for cdr3 hotspot in mutation
            if barCont[2].split('-')[0] == 'CDR3':
                if barCont[7] not in header:
                    header.append(barCont[7])
                isBarCode = list(barCode for barCode in mut_data if barCode['name'] == barCont[0].split('-')[0])
                if isBarCode == []:
                    mut_data.append({'name': barCont[0].split('-')[0], barCont[7]: 1})
                else:
                    if barCont[7] in isBarCode[0]:
                        isBarCode[0][barCont[7]] = isBarCode[0][barCont[7]] + 1
                    else:
                        isBarCode[0][barCont[7]] = 1
    
    # remove none from the data lists
    for i in range(len(mut_data)):
        if mut_data[i]['name'] == 'None':
            del mut_data[i]
            break 
    
    # write zeroes to barcodes without any data for particular mutations
    for mutation in header[1:]:
        for barCode in mut_data:
            if mutation not in barCode:
                barCode[mutation] = 0
    
    # write files to output directory
    wfName = os.path.join(outDir, "Mut_File_Based", "data", "mut_FileData.csv")
    with open(wfName, 'w', newline="") as writeFile:
        writer = csv.DictWriter(writeFile, fieldnames=header)
        writer.writeheader()
        for mutHS in mut_data:
            writer.writerow(mutHS)

def main():
    inFName = sys.argv[1]
    outDir = sys.argv[2]
    mut_CompareBarcodes(inFName, outDir)

main()

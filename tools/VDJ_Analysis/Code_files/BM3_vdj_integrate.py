# -*- coding: utf-8 -*-
"""
@author: scsac
"""

import csv

# integrate mutation data from imgt files to normalized matrix
def integrateMutations(inFileName, mutFileNames, outPath):
    fname = inFileName
    contigs = []
    with open(fname, 'r') as inFile:
        csvReader = csv.reader(inFile)
        contigs = next(csvReader)
    # read data from each mutation file for CDR1-3 and FR1-3    
    cdr3MutFName = mutFileNames[2]
    cdr3Content = []
    with open(cdr3MutFName, 'r') as cdr3File:
        cdr3Content = cdr3File.readlines()
    
    cdr3TotalMuts = [0]*len(contigs)
    cdr3HSMuts = [0]*len(contigs)
    cdr3NbHS = [0]*len(contigs)
    
    for cell in cdr3Content:
        cellData = cell.rstrip().split(",")
        cellName = cellData[0][:-3] + "_3"
        if cellName in contigs:
            cdr3TotalMuts[contigs.index(cellName)] = cellData[1]
            cdr3HSMuts[contigs.index(cellName)] = cellData[2]
            cdr3NbHS[contigs.index(cellName)] = cellData[3]
    
    cdr3TotalMuts[0] = "cdr3_total_muts"
    cdr3HSMuts[0] = "cdr3_hs_muts"
    cdr3NbHS[0] = "cdr3_nb_hs"
    
    cdr2MutFName = mutFileNames[1]
    cdr2Content = []
    with open(cdr2MutFName, 'r') as cdr2File:
        cdr2Content = cdr2File.readlines()
    
    cdr2TotalMuts = [0]*len(contigs)
    cdr2HSMuts = [0]*len(contigs)
    cdr2NbHS = [0]*len(contigs)
    
    for cell in cdr2Content:
        cellData = cell.rstrip().split(",")
        cellName = cellData[0][:-3] + "_3"
        if cellName in contigs:
            cdr2TotalMuts[contigs.index(cellName)] = cellData[1]
            cdr2HSMuts[contigs.index(cellName)] = cellData[2]
            cdr2NbHS[contigs.index(cellName)] = cellData[3]
    
    cdr2TotalMuts[0] = "cdr2_total_muts"
    cdr2HSMuts[0] = "cdr2_hs_muts"
    cdr2NbHS[0] = "cdr2_nb_hs"
    
    cdr1MutFName = mutFileNames[0]
    cdr1Content = []
    with open(cdr1MutFName, 'r') as cdr1File:
        cdr1Content = cdr1File.readlines()
    
    cdr1TotalMuts = [0]*len(contigs)
    cdr1HSMuts = [0]*len(contigs)
    cdr1NbHS = [0]*len(contigs)
    
    for cell in cdr1Content:
        cellData = cell.rstrip().split(",")
        cellName = cellData[0][:-3] + "_3"
        if cellName in contigs:
            cdr1TotalMuts[contigs.index(cellName)] = cellData[1]
            cdr1HSMuts[contigs.index(cellName)] = cellData[2]
            cdr1NbHS[contigs.index(cellName)] = cellData[3]
    
    cdr1TotalMuts[0] = "cdr1_total_muts"
    cdr1HSMuts[0] = "cdr1_hs_muts"
    cdr1NbHS[0] = "cdr1_nb_hs"
    
    fr3MutFName = mutFileNames[5]
    fr3Content = []
    with open(fr3MutFName, 'r') as fr3File:
        fr3Content = fr3File.readlines()
    
    fr3TotalMuts = [0]*len(contigs)
    fr3HSMuts = [0]*len(contigs)
    fr3NbHS = [0]*len(contigs)
    
    for cell in fr3Content:
        cellData = cell.rstrip().split(",")
        cellName = cellData[0][:-3] + "_3"
        if cellName in contigs:
            fr3TotalMuts[contigs.index(cellName)] = cellData[1]
            fr3HSMuts[contigs.index(cellName)] = cellData[2]
            fr3NbHS[contigs.index(cellName)] = cellData[3]
    
    fr3TotalMuts[0] = "fr3_total_muts"
    fr3HSMuts[0] = "fr3_hs_muts"
    fr3NbHS[0] = "fr3_nb_hs"
    
    fr2MutFName = mutFileNames[4]
    fr2Content = []
    with open(fr2MutFName, 'r') as fr2File:
        fr2Content = fr2File.readlines()
    
    fr2TotalMuts = [0]*len(contigs)
    fr2HSMuts = [0]*len(contigs)
    fr2NbHS = [0]*len(contigs)
    
    for cell in fr2Content:
        cellData = cell.rstrip().split(",")
        cellName = cellData[0][:-3] + "_3"
        if cellName in contigs:
            fr2TotalMuts[contigs.index(cellName)] = cellData[1]
            fr2HSMuts[contigs.index(cellName)] = cellData[2]
            fr2NbHS[contigs.index(cellName)] = cellData[3]
    
    fr2TotalMuts[0] = "fr2_total_muts"
    fr2HSMuts[0] = "fr2_hs_muts"
    fr2NbHS[0] = "fr2_nb_hs"
    
    fr1MutFName = mutFileNames[3]
    fr1Content = []
    with open(fr1MutFName, 'r') as fr1File:
        fr1Content = fr1File.readlines()
    
    fr1TotalMuts = [0]*len(contigs)
    fr1HSMuts = [0]*len(contigs)
    fr1NbHS = [0]*len(contigs)
    
    for cell in fr1Content:
        cellData = cell.rstrip().split(",")
        cellName = cellData[0][:-3] + "_3"
        if cellName in contigs:
            fr1TotalMuts[contigs.index(cellName)] = cellData[1]
            fr1HSMuts[contigs.index(cellName)] = cellData[2]
            fr1NbHS[contigs.index(cellName)] = cellData[3]
    
    fr1TotalMuts[0] = "fr1_total_muts"
    fr1HSMuts[0] = "fr1_hs_muts"
    fr1NbHS[0] = "fr1_nb_hs"
    # append data to the normalized matrix copied to the output directory
    wFName = outPath
    with open(wFName, 'a+', newline = "") as writeFile:
        csvWriter = csv.writer(writeFile)
        csvWriter.writerow(cdr3TotalMuts)
        csvWriter.writerow(cdr3HSMuts)
        csvWriter.writerow(cdr3NbHS)
        csvWriter.writerow(cdr2TotalMuts)
        csvWriter.writerow(cdr2HSMuts)
        csvWriter.writerow(cdr2NbHS)
        csvWriter.writerow(cdr1TotalMuts)
        csvWriter.writerow(cdr1HSMuts)
        csvWriter.writerow(cdr1NbHS)
        csvWriter.writerow(fr3TotalMuts)
        csvWriter.writerow(fr3HSMuts)
        csvWriter.writerow(fr3NbHS)
        csvWriter.writerow(fr2TotalMuts)
        csvWriter.writerow(fr2HSMuts)
        csvWriter.writerow(fr2NbHS)
        csvWriter.writerow(fr1TotalMuts)
        csvWriter.writerow(fr1HSMuts)
        csvWriter.writerow(fr1NbHS)
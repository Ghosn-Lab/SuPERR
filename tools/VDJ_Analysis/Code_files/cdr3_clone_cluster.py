# -*- coding: utf-8 -*-
"""
@author: scsac
"""

import sys

# reads the contig files and extracts the chain details
def readContigs(fileName):
    fileCont = []
    with open(fileName, 'r') as fileR:
        fileCont = fileR.readlines()
    
    cell_list = set()
    cdr3_list = {}
    IGH = ""
    IGKL = ""
    barC = ""
    for contig in fileCont[1:]:
        contDet = contig.rstrip().split(',')
        cell_list.add(contDet[0])
        if contDet[11] == 'True':
            if barC != contDet[0]:
                if barC != "":
                    if IGH != "" and IGKL != "":
                        if (IGH, IGKL) in cdr3_list:
                            cdr3_list[(IGH, IGKL)] += 1
                        else:
                            cdr3_list[(IGH, IGKL)] = 1
                    elif IGH != "":
                        if IGH in cdr3_list:
                            cdr3_list[IGH] += 1
                        else:
                            cdr3_list[IGH] = 1
                    else:
                        if IGKL in cdr3_list:
                            cdr3_list[IGKL] += 1
                        else:
                            cdr3_list[IGKL] = 1
                if contDet[5] == 'IGH':
                    barC = contDet[0]
                    IGH = contDet[12]
                    IGKL = ""
                elif contDet[5] == 'IGL' or contDet[5] == 'IGK':
                    barC = contDet[0]
                    IGKL = contDet[12]
                    IGH = ""
            else:
                if contDet[5] == 'IGH':
                    IGH = contDet[12]
                elif contDet[5] == 'IGL' or contDet[5] == 'IGK':
                    IGKL = contDet[12]
    if IGH != "" and IGKL != "":
        if (IGH, IGKL) in cdr3_list:
            cdr3_list[(IGH, IGKL)] += 1
        else:
            cdr3_list[(IGH, IGKL)] = 1
    elif IGH != "":
        if IGH in cdr3_list:
            cdr3_list[IGH] += 1
        else:
            cdr3_list[IGH] = 1
    else:
        if IGKL in cdr3_list:
            cdr3_list[IGKL] += 1
        else:
            cdr3_list[IGKL] = 1
    
    sort_cdr3 = sorted(cdr3_list.items(), key=lambda x: x[1], reverse=True)
    return sort_cdr3, len(cell_list)


def main():
    fileName = sys.argv[1] # get the input file from command line
    print(readContigs(fileName))
    
main()
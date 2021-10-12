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
    
    cdr3_list = {'IGH': set(), 'IGL': set(), 'IGK': set()}
    for contig in fileCont[1:]:
        barCode = contig.rstrip().split(',')
        if barCode[11] == 'True':
            cdr3_list[barCode[5]].add(barCode[12])
    
    return cdr3_list


def main():
    files_data = []
    # get the input file from command line
    for file in sys.argv[1:]:
        files_data.append(readContigs(file))
    
    IGH_Common = []
    IGL_Common = []
    IGK_Common = []
    
    # find intersection for each set with the next set
    i = 0
    while i < len(files_data)-1:
        IGH_Common = files_data[i]['IGH']
        IGL_Common = files_data[i]['IGL']
        IGK_Common = files_data[i]['IGK']
        
        IGH_Common = IGH_Common.intersection(files_data[i+1]['IGH'])
        IGL_Common = IGL_Common.intersection(files_data[i+1]['IGL'])
        IGK_Common = IGK_Common.intersection(files_data[i+1]['IGK'])
        
        i = i+1
    
    print(IGH_Common, IGK_Common, IGL_Common)
    
main()
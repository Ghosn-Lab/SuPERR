# -*- coding: utf-8 -*-
"""
@author: scsac
"""

import json, os
# generates the cdr3 aa treemap data
def cdr3AAData(inFileName, outDir):
    fname = inFileName
    fileCont = []
    with open(fname, 'r') as contig_file:
        fileCont = contig_file.readlines()
    # root level is the chain, check for amino acid presence, v, d and j presence and barcode addition
    cdr3Info = {'name': 'chain', 'children': []}
    for barcode in fileCont:
        barCont = barcode.split(",")
        if barCont[11] == 'True':
            isChain = list(chainN for chainN in cdr3Info['children'] if chainN['name'] == barCont[5])
            if isChain == []:
                cdr3Info['children'].append({'name': barCont[5], 'children': [{'name': barCont[12], 'value': 1,
                                                                               'barcodes': [barCont[0].split("-")[0]]}], 
                                             'V': [{'name': barCont[6], 'value': 1}], 'D': [{'name': barCont[7], 'value': 1}],
                                             'J': [{'name': barCont[8], 'value': 1}]})
            else:
                isAA = list(aaN for aaN in isChain[0]['children'] if aaN['name'] == barCont[12])
                isV = list(V for V in isChain[0]['V'] if V['name'] == barCont[6])
                if isV == []:
                    isChain[0]['V'].append({'name': barCont[6], 'value': 1})
                else:
                    isV[0]['value'] = isV[0]['value'] + 1
                isD = list(D for D in isChain[0]['D'] if D['name'] == barCont[7])
                if isD == []:
                    isChain[0]['D'].append({'name': barCont[7], 'value': 1})
                else:
                    isD[0]['value'] = isD[0]['value'] + 1
                isJ = list(J for J in isChain[0]['J'] if J['name'] == barCont[8])
                if isJ == []:
                    isChain[0]['J'].append({'name': barCont[8], 'value': 1})
                else:
                    isJ[0]['value'] = isJ[0]['value'] + 1
                if isAA == []:
                    isChain[0]['children'].append({'name': barCont[12], 'value': 1, 'barcodes': [barCont[0].split("-")[0]]})
                else:
                    isAA[0]['value'] = isAA[0]['value'] + 1
                    isAA[0]['barcodes'].append(barCont[0].split("-")[0])
    # write file to output directory
    wfName = os.path.join(outDir, "CDR3_AA_Treemap", "data", "cdr3_aa_data.json")
    with open(wfName, 'w') as wFile:
        wFile.write(json.dumps(cdr3Info))
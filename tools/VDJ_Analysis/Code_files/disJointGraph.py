# -*- coding: utf-8 -*-
"""
@author: scsac
"""

import os, sys, json

def disjointGraph(vdjContent):
    fileDict = {"nodes":[], "links":[]}
    for barCode in vdjContent:
        barCont = barCode.split(',')
        if barCont[11] == 'TRUE':
            # check for v gene in graph
            isV = list(vGene for vGene in fileDict['nodes'] if vGene['id'] == barCont[6])
            if isV == []:
                fileDict['nodes'].append({'id': barCont[6], 'group': "V", 'value': 1})
            else:
                isV[0]['value'] += 1
            # check for d gene in graph
            if barCont[7] != 'None':
                isD = list(dGene for dGene in fileDict['nodes'] if dGene['id'] == barCont[7])
                if isD == []:
                    fileDict['nodes'].append({'id': barCont[7], 'group': "D", 'value': 1})
                else:
                    isD[0]['value'] += 1
                # check for v-d link in the graph
                isDLink = list(dLink for dLink in fileDict['links'] if dLink['source'] == barCont[6] and dLink['target'] == barCont[7])
                if isDLink == []:
                    fileDict['links'].append({'source': barCont[6], 'target': barCont[7], 'value': 1})
                else:
                    isDLink[0]['value'] += 1
            # check for j gene in graph
            isJ = list(jGene for jGene in fileDict['nodes'] if jGene['id'] == barCont[8])
            if isJ == []:
                fileDict['nodes'].append({'id': barCont[8], 'group': "J", 'value': 1})
            else:
                isJ[0]['value'] += 1
            # check for v-j link in the graph
            isJLink = list(jLink for jLink in fileDict['links'] if jLink['source'] == barCont[6] and jLink['target'] == barCont[8])
            if isJLink == []:
                fileDict['links'].append({'source': barCont[6], 'target': barCont[8], 'value': 1})
            else:
                isJLink[0]['value'] += 1
    return fileDict

def main():
    fileContent = []
    fName = sys.argv[1]
    with open(fName, 'r') as fRead:
        fileContent = fRead.readlines()
    writeContent = disjointGraph(fileContent[1:])
    wFName = os.path.join(sys.argv[2], 'vdjRecombination.json')
    with open(wFName, 'w') as fWrite:
        fWrite.write(json.dumps(writeContent))

main()
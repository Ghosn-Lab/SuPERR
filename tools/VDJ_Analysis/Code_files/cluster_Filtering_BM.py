# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 14:50:26 2020

@author: scsac
"""

import sys, os

# parse the cluster files and returns the lists for CD138 pos and neg clusters for BM2 and BM3 samples
def parseClusterFile(fileName):
    clusterFileCont = []
    with open(fileName, 'r') as fileReader:
        clusterFileCont = fileReader.readlines()
    
    CD138_BM2_0, CD138_BM2_1, CD138_BM3_0, CD138_BM3_1 = [], [], [], []
    for contig in clusterFileCont[1:]:
        contigDet = contig.rstrip().split(',')
        contigBC = contigDet[0].split('_')
        if contigDet[1] == "0":
            if contigBC[1] == "BM2":
                CD138_BM2_0.append(contigBC[0])
            else:
                CD138_BM3_0.append(contigBC[0])
        else:
            if contigBC[1] == "BM2":
                CD138_BM2_1.append(contigBC[0])
            else:
                CD138_BM3_1.append(contigBC[0])

    return CD138_BM2_0, CD138_BM2_1, CD138_BM3_0, CD138_BM3_1

# generates the cluster based CD138 cell matrices for BM2 and BM3
def main():
    clusterFile1 = sys.argv[1]
    clusterFile2 = sys.argv[2]
    contig_file1 = sys.argv[3]
    contig_file2 = sys.argv[4]
    outDir = sys.argv[5]
    
    CD138_pos_BM2_0, CD138_pos_BM2_1, CD138_pos_BM3_0, CD138_pos_BM3_1 = parseClusterFile(clusterFile1)
    CD138_neg_BM2_0, CD138_neg_BM2_1, CD138_neg_BM3_0, CD138_neg_BM3_1 = parseClusterFile(clusterFile2)
    
    contig_fileCont1 = []
    contig_fileCont2 = []
    with open(contig_file1, 'r') as fileR1:
        contig_fileCont1 = fileR1.readlines()
    with open(contig_file2, 'r') as fileR2:
        contig_fileCont2 = fileR2.readlines()
    
    # append the headers for the contig files
    CD138_pos_cluster_0 = [contig_fileCont1[0]]
    CD138_pos_cluster_1 = [contig_fileCont1[0]]
    CD138_neg_cluster_0 = [contig_fileCont1[0]]
    CD138_neg_cluster_1 = [contig_fileCont1[0]]
    
    # filter the BM2 contig file by cluster and CD138
    for contigs in contig_fileCont1[1:]:
        barCode = contigs.split(',')[0].split('-')[0]
        if barCode in CD138_pos_BM2_0:
            CD138_pos_cluster_0.append(contigs)
        elif barCode in CD138_neg_BM2_0:
            CD138_neg_cluster_0.append(contigs)
        elif barCode in CD138_pos_BM2_1:
            CD138_pos_cluster_1.append(contigs)
        elif barCode in CD138_neg_BM2_1:
            CD138_neg_cluster_1.append(contigs)
        else:
            continue

    # filter the BM3 contig file by cluster and CD138
    for contigs in contig_fileCont2[1:]:
        barCode = contigs.split(',')[0].split('-')[0]
        if barCode in CD138_pos_BM3_0:
            CD138_pos_cluster_0.append(contigs)
        elif barCode in CD138_neg_BM3_0:
            CD138_neg_cluster_0.append(contigs)
        elif barCode in CD138_pos_BM3_1:
            CD138_pos_cluster_1.append(contigs)
        elif barCode in CD138_neg_BM3_1:
            CD138_neg_cluster_1.append(contigs)
        else:
            continue
    
    # write the cluster and CD138 filtered sets to files
    outFile1 = os.path.join(outDir, "CD138_pos_cluster0.csv")
    outFile2 = os.path.join(outDir, "CD138_pos_cluster1.csv")
    outFile3 = os.path.join(outDir, "CD138_neg_cluster0.csv")
    outFile4 = os.path.join(outDir, "CD138_neg_cluster1.csv")
    with open(outFile1, 'w') as wFile:
        wFile.writelines(CD138_pos_cluster_0)
    with open(outFile2, 'w') as wFile:
        wFile.writelines(CD138_pos_cluster_1)
    with open(outFile3, 'w') as wFile:
        wFile.writelines(CD138_neg_cluster_0)
    with open(outFile4, 'w') as wFile:
        wFile.writelines(CD138_neg_cluster_1)

main()
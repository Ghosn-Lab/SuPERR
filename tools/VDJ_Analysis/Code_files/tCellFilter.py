# -*- coding: utf-8 -*-
"""
@author: scsac
"""

import csv

def tCellFilter(fName, outFName):
    fileContents = []
    with open(fName, 'r') as fileReader:
        fileContents = fileReader.readlines()
    
    header = ['barcode', 'length_A', 'length_B', 'v_gene_A', 'j_gene_A', 'c_gene_A', 'v_gene_B', 'j_gene_B', 'c_gene_B',
              'cdr3_A', 'cdr3_B', 'cdr3_nt_A', 'cdr3_nt_B', 'reads_A', 'reads_B', 'umis_A', 'umis_B', 'clonotype_id']
    
    filtered_contigs = []
    for barcode in fileContents[1:]:
        contig = barcode.rstrip().split(',')
        if contig[11] == 'True':
            isBarcode = list(contigs for contigs in filtered_contigs if contigs['barcode'] == contig[0])
            if isBarcode == []:
                if contig[5] == 'TRA':
                    filtered_contigs.append({'barcode': contig[0], 'length_A': contig[4], 'length_B': 'None', 
                                             'v_gene_A': contig[6], 'j_gene_A': contig[8], 'c_gene_A': contig[9],
                                             'v_gene_B': 'None', 'j_gene_B': 'None', 'c_gene_B': 'None',
                                             'cdr3_A': contig[12], 'cdr3_B': 'None', 'cdr3_nt_A': contig[13],
                                             'cdr3_nt_B': 'None', 'reads_A': contig[14], 'reads_B': 'None',
                                             'umis_A': contig[15], 'umis_B': 'None', 'clonotype_id': contig[16]})
                else:
                    filtered_contigs.append({'barcode': contig[0], 'length_A': 'None', 'length_B': contig[4], 
                                             'v_gene_A': 'None', 'j_gene_A': 'None', 'c_gene_A': 'None',
                                             'v_gene_B': contig[6], 'j_gene_B': contig[8], 'c_gene_B': contig[9],
                                             'cdr3_A': 'None', 'cdr3_B': contig[12], 'cdr3_nt_A': 'None',
                                             'cdr3_nt_B': contig[13], 'reads_A': 'None', 'reads_B': contig[14],
                                             'umis_A': 'None', 'umis_B': contig[15], 'clonotype_id': contig[16]})
            else:
                if contig[5] == 'TRA':
                    isBarcode[0]['length_A']= contig[4]
                    isBarcode[0]['v_gene_A']= contig[6]
                    isBarcode[0]['j_gene_A']= contig[8]
                    isBarcode[0]['c_gene_A']= contig[9]
                    isBarcode[0]['cdr3_A']= contig[12]
                    isBarcode[0]['cdr3_nt_A']= contig[13]
                    isBarcode[0]['reads_A']= contig[14]
                    isBarcode[0]['umis_A']= contig[15]
                else:
                    isBarcode[0]['length_B']= contig[4]
                    isBarcode[0]['v_gene_B']= contig[6]
                    isBarcode[0]['j_gene_B']= contig[8]
                    isBarcode[0]['c_gene_B']= contig[9]
                    isBarcode[0]['cdr3_B']= contig[12]
                    isBarcode[0]['cdr3_nt_B']= contig[13]
                    isBarcode[0]['reads_B']= contig[14]
                    isBarcode[0]['umis_B']= contig[15]
        
    with open(outFName, 'w', newline="") as fileWriter:
        writer = csv.DictWriter(fileWriter, fieldnames=header)
        writer.writeheader()
        for contig in filtered_contigs:
            writer.writerow(contig)

fName = 'C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\VDJ\Danny\\p20090-s004_4_051920_DC_L2_5GEX_D1\\filtered_contig_annotations.csv'
outFName = 'C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\VDJ\Danny\\p20090-s004_4_051920_DC_L2_5GEX_D1\\results.csv'
tCellFilter(fName, outFName)
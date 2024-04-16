#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mabin
"""

####Bd-RPCE 
###methy gene select 
from Bio import SeqIO
import pandas as pd 
import os 
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-in',help="Input E. coli sequence [required]")
parser.add_argument('-identity',help="Blastn identity cutoff (default 80)")
parser.add_argument('-cover',help="Blastn coverage cutoff (default 80)")
parser.add_argument('-o',help="Output location [required]")
parser.add_argument('-word_size',help="Blastn searching word size (default 11)")
parser.add_argument('-threads',help="Blastn thread number (default 1)")
parser.parse_args()


####input paraments
if len(sys.argv) == 1:
    print('Please input the paraments')
    sys.exit(0)
####  
  
#intput location
try:
    index = sys.argv.index('-in')+1
    seq_location  = sys.argv[index]
    for seq_record in SeqIO.parse(seq_location, 'fasta'):
        seq_record.seq
    print('The input sequence is %s' %seq_location)
except:
    print('Please input correct input sequence')
    sys.exit(0)
    
#output location
try:
    index = sys.argv.index('-o')+1
    output_location  = sys.argv[index]
    if os.path.isdir(output_location) == False:
        1/0
    print('The output location is %s' %output_location)
except:
    print('Please input correct output directory')
    sys.exit(0)
    
#identity
try:
    index = sys.argv.index('-identity')+1
    identity_cutoff  = float(sys.argv[index])
    print('Blastn identity cutoff is %s percent' % str(identity_cutoff))
except:
    identity_cutoff = 80
    print('Blastn identity cutoff is 80 percent')

#cover
try:
    index = sys.argv.index('-cover')+1
    cover_cutoff  = float(sys.argv[index])
    print('Blastn cover cutoff is %s percent' % str(cover_cutoff))
except:
    cover_cutoff = 80
    print('Blastn cover cutoff is 80 percent')

#word size
try:
    index = sys.argv.index('-word_size')+1
    word_size = int(sys.argv[index])
    print('Blastn word size is %s' % str(word_size))
except:
    word_size = 11
    print('Blastn word size is 11 percent')

#cover
try:
    index = sys.argv.index('-threads')+1
    threads  = int(sys.argv[index])
    print('Blastn threads is %s' % str(threads))
except:
    threads = 1
    print('Blastn threads is 1')

#############

total_seq_id = []
total_seq = []
for seq_record in SeqIO.parse(seq_location,'fasta'):
    total_seq_id.append(seq_record.id)    
    total_seq.append(seq_record.seq)
    
methy_gene = []
methy_gene_length = []
for seq_record in SeqIO.parse('./methy_ref.fasta', 'fasta'):
    methy_gene.append(seq_record.id)
    methy_gene_length.append(len(seq_record.seq))

os.system('blastn -query %s -db ./methy_ref -outfmt 6 -out ./tmp_methy_gene -word_size %s -num_threads %s' % (seq_location,word_size,threads))
total_seq_info =[]
out_seq_info = []
for i in range(len(total_seq_id)):              
    blast_out = pd.read_table('./tmp_methy_gene',header=None)
    methy_gene_index = 0
    tmp_gene_seq = []
    blast_out = blast_out[blast_out[0]==total_seq_id[i]]
    for j in range(len(methy_gene)):
        blast_out_select = blast_out[blast_out[1]==methy_gene[j]]
        if len(blast_out_select) == 0:
            print(methy_gene[j])
            continue
        blast_out_select = blast_out_select.sort_values(by=11,ascending=False).reset_index(drop=True)
        identity = blast_out_select[2][0]
        cover = blast_out_select[3][0]/methy_gene_length[j] * 100
        # cover=0.8
        if identity>=identity_cutoff and cover>=cover_cutoff:
            start = blast_out_select[6][0]
            end = blast_out_select[7][0]
            gene_start = blast_out_select[8][0]
            gene_end = blast_out_select[9][0]
            chrom = blast_out_select[0][0]
            
            index = total_seq_id.index(chrom)
            strain_seq = total_seq[index]
  
            if gene_start < gene_end:
                gene_sequence = strain_seq[start:end]
            else:
                gene_sequence = strain_seq[start:end].reverse_complement()
            methy_gene_index+=1
            tmp_gene_seq.append([methy_gene[j],gene_sequence])
            
    if methy_gene_index != 20:
        out_seq_info.append([total_seq_id[i],methy_gene_index])
        # break
    else:
        tmp_sequence = ''
        for k in range(len(tmp_gene_seq)):
            tmp_sequence = tmp_sequence+ tmp_gene_seq[k][1]
        total_seq_info.append([total_seq_id[i],tmp_sequence])

output_result = output_location+'HMC_gene.fasta'      
with open(output_result,'w') as f:
    for i in range(len(total_seq_info)):
        f.write('>'+total_seq_info[i][0]+'\n'+str(total_seq_info[i][1])+'\n')
        
output_NA = output_location+'NA_strain.info'
if len(out_seq_info) == 0:
    pass
else:
    pd.DataFrame(out_seq_info).to_csv(output_NA,header=['ID','HMCgene number'],index=None)













#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
BdRPC_CNS.py

                    Bd-RPC 
(Bases dependent Rapid Phylogenetic Clustering)

            CLUSTERING NEW SEQUENCES
                                
                                Author: Ma Bin
'''
#####Clustering new sequences function
def calcuate_bases_frequency(aligned_seq_location):
    from Bio import SeqIO
    A = []
    C = []
    G = []
    T = []
    seq_len = 0 ##test whether aligned 
    seq_id = []
    sequence_out = []
    for seq_record in SeqIO.parse(aligned_seq_location,'fasta'):
        seq_id.append(seq_record.id)
        sequence = seq_record.seq.lower()  ####change (A C G T) to (a c g t)
        A.append(sequence.count('a'))
        C.append(sequence.count('c'))
        G.append(sequence.count('g'))
        T.append(sequence.count('t'))
        
        ###test aligned sequence
        if seq_len == 0:
            seq_len = len(sequence)
        
        if seq_len != len(sequence):
            exit(print('Please input aligned sequences'))
            
        sequence_out.append(sequence)
        ###########################
    freq_A = sum(A)/(sum(A)+sum(C)+sum(G)+sum(T))
    freq_C = sum(C)/(sum(A)+sum(C)+sum(G)+sum(T))
    freq_G = sum(G)/(sum(A)+sum(C)+sum(G)+sum(T))
    freq_T = sum(T)/(sum(A)+sum(C)+sum(G)+sum(T))
    
    print ('Frequency A : '+ str(freq_A) +'\n'+
           'Frequency C : '+ str(freq_C) +'\n'+
           'Frequency G : '+ str(freq_G) +'\n'+
           'Frequency T : '+ str(freq_T))
    return [freq_A,freq_C,freq_G,freq_T],seq_id,sequence_out

def bases_convert(pi, sequence, convert_rule_location = '' ):
    import numpy as np
    from Bio import SeqIO
    if convert_rule_location == '':
        A = np.array([1,0,0,0,1,0])* (1-pi[0])
        C = np.array([0,1,0,0,0,1])* (1-pi[1])
        G = np.array([0,0,1,0,1,0])* (1-pi[2])
        T = np.array([0,0,0,1,0,1])* (1-pi[3])
    else:
        convert_rule = np.loadtxt(convert_rule_location ,delimiter = ',',encoding = 'utf-8-sig') ###sort by A C G T
        A = convert_rule[0,:]
        C = convert_rule[1,:]
        G = convert_rule[2,:]
        T = convert_rule[3,:]
    
    R = (A + G)/2
    Y = (C + T)/2
    S = (G + C)/2
    W = (A + T)/2
    K = (G + T)/2
    M = (A + C)/2
    B = (C + G + T)/3
    D = (A + G + T)/3
    H = (A + C + T)/3
    V = (A + C + G)/3
    gap = N = (A + C + G + T)/4
    
    convert_dict = {'a':A,'c':C,'g':G,'t':T,'-':gap,'r':R,'y':Y,'s':S,'w':W,'k':K,'m':M,'b':B,'d':D,'h':H,'v':V,'n':N}
    seq_change_matrix = []
    
    for i in range(len(sequence)):
        tmp_seq = [convert_dict[k] if k in convert_dict else k for k in list(sequence[i])]
              
        tmp_seq = np.array(tmp_seq)
        tmp_seq = tmp_seq.reshape(1,tmp_seq.shape[0]*tmp_seq.shape[1])
    
        seq_change_matrix.append(tmp_seq[0])
        
    seq_change_matrix = np.array(seq_change_matrix)
    
    return seq_change_matrix,[A,C,G,T]
    
def PCA_improved(seq_change_matrix,PCA_components = 'max'):
    
    from sklearn.decomposition import PCA
    import numpy as np
    
    seq_change_matrix = np.array(seq_change_matrix)
    
    if PCA_components == 'max':
        PCA_components = seq_change_matrix.shape[0]
    else:
        PCA_components = int(PCA_components)
    
    
    pca = PCA(n_components=PCA_components)
    pca.fit(seq_change_matrix)
    seq_change_matrix_PCA  = pca.fit_transform(seq_change_matrix)
    
    #print ('PCA explained variance = ' + str(sum(pca.explained_variance_ratio_)))
    
    return seq_change_matrix_PCA
 
def select_tree_id(tree):
    tree_id = []
    total_tree_id = tree.get_terminals()
    for i in range(len(total_tree_id)):
        tree_id.append(total_tree_id[i].name)
    return tree_id
           
def add_aligned_sequence(add_seq_location, aligned_seq_location,output_location = '',thread = 1):
    import os 
    if output_location == '':
        output_location = aligned_seq_location + '.combine'
        
    os.system("mafft  --thread -1 --add %s --reorder %s > %s" %( add_seq_location, aligned_seq_location, output_location))
    
    return 

def gap_t_test(add_seq_location, aligned_seq_location, output_location = '',IDfold=1.1):
    from Bio import SeqIO
    import numpy as np
    new_seqid = []
    gap = 1
    
    ###add seq id
    add_seq_id = []
    for seq_record in SeqIO.parse(add_seq_location,'fasta'):
        add_seq_id.append(seq_record.id)
    
    ###change seq format -- aligned before
    before_id = []
    before_information = []
    for seq_record in SeqIO.parse(aligned_seq_location,'fasta'):
        tmp=[]
        before_id.append(seq_record.id)
        for i in range(len(seq_record.seq)):
            if seq_record.seq[i] == '-':
                tmp.append(gap)
            else:
                tmp.append(0)
        before_information.append(tmp)
        
    before_information = np.array(before_information)
    
    ###change seq format -- aligned after
    if output_location == '':
        output_location = aligned_seq_location + '.combine'
    
    after_id = []
    after_information = []
    for seq_record in SeqIO.parse(output_location,'fasta'):
        tmp=[]
        after_id.append(seq_record.id)
        for i in range(len(seq_record.seq)):
            if seq_record.seq[i] == '-':
                 tmp.append(gap)
            else:
                tmp.append(0)
        after_information.append(tmp)
        
    after_information = np.array(after_information)
    
    ###match
    before_match_information = []
    after_match_information = []
    for i in range(len(before_id)):
        before_match_information.append(before_information[i,:])
        after_match_information.append(after_information[after_id.index(before_id[i]),:])
    
    before_match_information = np.array(before_match_information)
    after_match_information = np.array(after_match_information)
    
    ###distinguish different sequence
    before_test = sum(before_match_information.T)
    out_seq_id = []
    median_list = []
    for i in range(len(add_seq_id)):
        
        add_seq_information = after_information[after_id.index(add_seq_id[i]),:]
        
        after_match_information_combine = np.vstack((after_match_information,add_seq_information))
        delete_same_gap = len(np.where(sum(after_match_information_combine)==after_match_information_combine.shape[0])[0])
        
        after_test = sum(after_match_information.T) - delete_same_gap
            
        if np.median(after_test) <= IDfold* np.median(before_test):
            new_seqid.append(add_seq_id[i])
            print (str(add_seq_id[i])+'\t'+'IN Database')
        else:
            print (str(add_seq_id[i])+'\t'+'OUT Database')
            print (np.median(after_test))
            print (np.median(before_test))
            out_seq_id.append(add_seq_id[i])
        
        median_list.append(np.median(after_test))
    return [new_seqid,out_seq_id,median_list]

def clustering_information_tree(database_location,aligned_seq_location,combine_seq_location,new_seqid,convert_rule_location,identity_cutoff=0.8,density_fold=1.5,pdistance=''):
    from Bio import SeqIO
    import numpy as np
    from scipy.spatial.distance import pdist, squareform

    ###############
    database = []
    with open(database_location,'r') as f:
        for i in f.readlines():
            i = i.split('\n')
            database.append(i[0])
    
    PCA_switch = database[0].split(',')[0].split('/')[0]
    PCA_components = database[0].split(',')[0].split('/')[1]
    distance_exponent = float(database[0].split(',')[1])
    
    ##remove low identity cluster
    for i in range(len(database)-1,0,-1):
        if float(database[i].split('/')[1]) < float(identity_cutoff):
            database.remove(database[i])
    
    ####change seqid sequence  order
    ##input combine seq_id sequence
    combine_seq_id = []
    combine_sequence = []
    
    for seq_record in SeqIO.parse(combine_seq_location,'fasta'):
        combine_seq_id.append(seq_record.id)
        combine_sequence.append(seq_record.seq)
    ##input database seq_id sequence
    seq_id = []
    sequence = []
    for seq_record in SeqIO.parse(aligned_seq_location, 'fasta'):
        seq_id.append(seq_record.id)
        sequence.append(seq_record.seq)
        
    ##resort seq_id sequence. database+new
    final_seq_id = [''] * len(seq_id)
    final_sequence = [''] * len(sequence)
    for i in range(len(combine_seq_id)):
        try:
            index = seq_id.index(combine_seq_id[i])
            final_seq_id[index] = combine_seq_id[i]
            final_sequence[index] = combine_sequence[i]
        except:
            pass
        try:
            index = new_seqid.index(combine_seq_id[i])
            final_seq_id.append(combine_seq_id[i])
            final_sequence.append(combine_sequence[i])
        except:
            pass
    #######################################################      
    
    ###convert sequence by database convert matrix
    if len(pdistance) != 0:
        distance_matrix=np.zeros([len(final_seq_id),len(final_seq_id)])
        for i in range(len(final_seq_id)):
            for j in range(i+1,len(final_seq_id)):
                index = list(pdistance['Unnamed: 0']).index(final_seq_id[j])
                distance_matrix[i,j]=pdistance[final_seq_id[i]][index]
                distance_matrix[j,i]=pdistance[final_seq_id[i]][index]
    else:
        seq_change_matrix,convert_matrix= bases_convert(0,sequence = final_sequence,convert_rule_location=convert_rule_location)
        
        ###PCA increasing
        ##PCA_switch PCA_components
        if PCA_switch == 'PCA_on':
            seq_change_matrix = PCA_improved(seq_change_matrix,PCA_components)
        else:
            pass
        
        ###calcuate distance matrix
        ##distance_exponent
        if  distance_exponent == 2:
            distance_matrix = pdist(seq_change_matrix,'euclidean')
            distance_matrix = squareform(distance_matrix)
        elif distance_exponent == 1:
            distance_matrix = pdist(seq_change_matrix,'cityblock')
            distance_matrix = squareform(distance_matrix)
        else:
            distance_matrix = pdist(seq_change_matrix,'minkowski',p=distance_exponent)
            distance_matrix = squareform(distance_matrix)
            
    ####clustering new sequence
    database_cutoff = len(seq_id)
    clustering_output = []
    clustering_density = []
    clustering_id = []
    for i in range (database_cutoff,len(final_seq_id)):
        seach_distance_matrix  = distance_matrix[i][:database_cutoff]
        min_location =np.where(seach_distance_matrix == np.min(seach_distance_matrix))[0][0] 
        clustering_output.append(min_location)
        clustering_density.append(seach_distance_matrix[min_location])
        clustering_id.append(final_seq_id[i])
        
    ####search in database 
    clustering_DBsearch = []
    clustering_DBsearch_id = []
    clustering_DBsearch_total = []
    for i in range(len(clustering_output)):
        tmp_search_out = []
        #tmp_search_out_total=[]
        for j in range(1,len(database)):
            try:
                database[j].split('/')[-1].split(',').index(str(clustering_output[i]))
                #tmp_search_out_total.append(database[j].split('/')[0])
                if clustering_density[i] < density_fold * float(database[j].split('/')[2]):
                    tmp_search_out.append(database[j].split('/')[0])
            except:
                pass
        clustering_DBsearch_id.append(clustering_id[i])
        if tmp_search_out == []:
            tmp_search_out = ['empty']
        clustering_DBsearch.append(tmp_search_out)
        #clustering_DBsearch_total.append(tmp_search_out_total)
    ######select min-tree 
    clustering_result = []
    for i in range(len(clustering_DBsearch)):
        clustering_result.append([clustering_DBsearch_id[i],clustering_DBsearch[i][-1]])
    
    clustering_result_total = []
    # for i in range(len(clustering_DBsearch)):
    #     clustering_result_total.append([clustering_DBsearch_id[i],clustering_DBsearch_total[i][:]])
    
    return [clustering_result,clustering_result_total,clustering_density,clustering_output]

def add_tree(clustering_result,ML_tree_location,aligned_seq_location,combine_seq_location,threads=1):
    from Bio import SeqIO
    from Bio import Phylo
    from Bio.Phylo.BaseTree import Tree
    import pandas as pd
    import os 
    from Bio.SeqRecord import SeqRecord
    from io import StringIO
    import random
    import numpy as np
    ###def function
    def read_tree(tree,tree_location):
        #tree = tree.clade
        tree_location = list(tree_location)
        for i in range (len(tree_location)):
            tree = tree.clades[int(tree_location[i])]
        return tree
    
    def select_tree_id(tree):
        tree_id = []
        total_tree_id = tree.get_terminals()
        for i in range(len(total_tree_id)):
            tree_id.append(total_tree_id[i].name)
        return tree_id
    
    tree = Phylo.read(ML_tree_location, 'newick')
    tree.root_at_midpoint()
    tree= tree.clade
    
    old_tree = Phylo.read(ML_tree_location, 'newick')
    old_tree.root_at_midpoint()
    old_tree= old_tree.clade
    
    ##input combine seq_id sequence
    combine_seq_id = []
    combine_sequence = []
    
    for seq_record in SeqIO.parse(combine_seq_location,'fasta'):
        combine_seq_id.append(seq_record.id)
        combine_sequence.append(seq_record.seq)
    
    ###remove same location
    location_list_index = list(pd.value_counts(pd.DataFrame(clustering_result)[1]).index)
    clustering_result_pd = pd.DataFrame(clustering_result)

    add_seq_id_list = []
    for i in range(len(location_list_index)):
        add_seq_id_list.append(list(clustering_result_pd[clustering_result_pd[1]==location_list_index[i]][0]))
    
    for j in reversed(range(len(location_list_index))):
        for i in range(len(location_list_index)):
            if i == j:
                continue
            try:
                if location_list_index[j].find(location_list_index[i]) == 0:
                    
                    for k in range(len(add_seq_id_list[j])):
                        add_seq_id_list[i].append(add_seq_id_list[j][k])
                    
                    del add_seq_id_list[j]
                    del location_list_index[j]
            except:
                pass
    #####################
    time_info = []
    ##merge same tree location 
    for i in range(len(location_list_index)):
        tree_location =location_list_index[i]
        print(tree_location)
        add_seq_id = add_seq_id_list[i]
        db_seq_id = select_tree_id(read_tree(old_tree,tree_location))
        final_add_seq_id = add_seq_id + db_seq_id
        ##write fasta to make tree
        tmp_out_fasta = []
        for j in range(len(combine_seq_id)):
            if combine_seq_id[j] in final_add_seq_id:
                tmp_out_fasta.append(SeqRecord(combine_sequence[j],id=combine_seq_id[j],description=''))
    
        
        if len(tmp_out_fasta) == 2:
            new_tree = Phylo.read(StringIO("(%s,%s)" %(tmp_out_fasta[0].id,tmp_out_fasta[1].id)), "newick")
        elif len(tmp_out_fasta) == 3:
            SeqIO.write(tmp_out_fasta, './tmp_out.fasta', 'fasta')
            os.system('iqtree -s ./tmp_out.fasta -T %s ' %threads)
            new_tree = Phylo.read('./tmp_out.fasta.treefile', 'newick')
            new_tree.root_at_midpoint()
            os.system('rm ./tmp_out.*')
        else:
            SeqIO.write(tmp_out_fasta, './tmp_out.fasta', 'fasta')
            os.system('iqtree -s ./tmp_out.fasta -B 1000 -T %s ' % threads)
            new_tree = Phylo.read('./tmp_out.fasta.treefile', 'newick')
            new_tree.root_at_midpoint()
            os.system('rm ./tmp_out.*')
        
        # Phylo.write(new_tree.clade,'./%s' % tree_location,'newick')
        
        # new_location = './%s' % tree_location+'.csv'
        # with open(new_location ,'w') as f:
        #     for h in range(len(add_seq_id)):
        #         f.write(str(add_seq_id[h])+'\n')
            
        tree_location_change = 'tree'
        for j in range(len(list(tree_location))):
            tree_location_change += '.clades['+str(list(tree_location)[j])+']'
        
        exec(tree_location_change + ' = new_tree.clade')
        exec(tree_location_change + '.branch_length= 0')
        exec(tree_location_change + '.confidence= 100')
        
        ##calcuate branch length
        if list(location_list_index[i])[-1] =='0':
            other_location = ''.join(list(location_list_index[i])[0:-1]) + '1'
        else:
            other_location = ''.join(list(location_list_index[i])[0:-1]) + '0'
        
        ###select other id 
        if len(select_tree_id(read_tree(old_tree,other_location))) > 5:
            other_id = random.sample(select_tree_id(read_tree(old_tree,other_location)),5)
        else:
            other_id = select_tree_id(read_tree(old_tree,other_location))
        ###select location id 
        if len(select_tree_id(read_tree(old_tree,location_list_index[i]))) > 5:
            in_id = random.sample(select_tree_id(read_tree(old_tree,location_list_index[i])),5)
        else:
            in_id = select_tree_id(read_tree(old_tree,location_list_index[i]))
        
        branch = []
        for h in range(len(other_id)):
            for k in range(len(in_id)):
                branch.append(old_tree.distance(other_id[h],in_id[k]) - tree.distance(other_id[h],in_id[k]))
        
        branch_length = np.mean(branch)
        
        #print(branch_length)
        
        exec(tree_location_change + '.branch_length= %s' % branch_length)
    return tree
        

import sys
import os
from Bio import SeqIO
from Bio import Phylo
import numpy as np
import pandas as pd
import argparse

####parament set
parser = argparse.ArgumentParser()
parser.add_argument('-align',help="Location of aligned sequences. (required) [no punctuation mark: '/' or ',']")
parser.add_argument('-db',help="Location of Bd-RPC database. (required)")
parser.add_argument('-new',help="Location of new sequences. (required)")
parser.add_argument('-o',help="Directory to store the result. (required)")
parser.add_argument('-IDfold',help="The fold of median value in Indel Test. (default: 1.1)")
parser.add_argument('-phy_information',help="Location of phylogentics tree. (if the tree is available, the new sequences will be inserted into the phylogenetic tree)")
parser.add_argument('-identity_cutoff',help="The cutoff value of clusters' identity. (0~1, default: 0.8)")
parser.add_argument('-density_fold',help="The fold of clusters' density for new samples clustering. (default: 1.5)")
parser.add_argument('-threads',help="Threads of mafft align and iqtree. (int, default: 1)")
parser.parse_args()

####input paraments
if len(sys.argv) == 1:
    print('Please input the paraments')
    sys.exit(0)
##Check paraments
#aligned sequences
seq_number = 0
seq_id = []
try:
    index = sys.argv.index('-align')+1
    aligned_seq_location = sys.argv[index]
    for seq_record in SeqIO.parse(aligned_seq_location, 'fasta'):
        seq_length = len(seq_record.seq)
        break
    for seq_record in SeqIO.parse(aligned_seq_location, 'fasta'):
        seq_id.append(seq_record.id)
        seq_number+=1
        if len(seq_record.seq) != seq_length:
            1/0
    seq_id.sort()
except:
    print('Please input correct aligned sequences')
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

#new sequences
try:
    index = sys.argv.index('-new')+1
    addition_seq_location  = sys.argv[index]
    new_seqid = []
    for seq_record in SeqIO.parse(addition_seq_location, 'fasta'):
        new_seqid.append(seq_record.id)
    if new_seqid == []:
        1/0
except:
    print('Please input correct new sequences location (fasta format)')
    sys.exit(0)
    
#combine sequences
try:
    index = sys.argv.index('-combine')+1
    combine_location  = sys.argv[index]
    combine_seqid = []
    for seq_record in SeqIO.parse(combine_location, 'fasta'):
        combine_seqid.append(seq_record.id)
    if len(combine_seqid) != len(seq_id)+len(new_seqid):
        combine_location = ''
        print('Combine aligned sequences are not correct. The Mafft will be used to aligned the sequences.')
except:
    combine_location = ''
    print('Combine aligned sequences are not correct. The Mafft will be used to aligned the sequences.')
#database location
try:
    index = sys.argv.index('-db')+1
    database_location  = sys.argv[index]
    database = []
    with open(database_location+'.match','r') as f:
        for i in f.readlines():
            i = i.split('\n')
            database.append(i[0])
    
    PCA_switch = database[0].split('/')[0]
    PCA_components = database[0].split(',')[0].split('/')[1]
    distance_exponent = float(database[0].split(',')[1])
    
    # convert_rule = np.loadtxt(database_location+'.convert' ,delimiter = ',',encoding = 'utf-8-sig') ###sort by A C G T
except:
    print('Please input correct database location (ps: /../../database)')
    sys.exit(0)

#IDfold
try:
    index = sys.argv.index('-IDfold')+1
    IDfold  = sys.argv[index]
except:
    IDfold = 1.1


try:
    IDfold = float(IDfold)
except:
    print('Please input correct IDfold (float, default: 1.1)')
    sys.exit(0)


#phy_information
try:
    index = sys.argv.index('-phy_information')+1
    phy_information  = sys.argv[index]
except:
    phy_information = ''

try:
    if phy_information == '':
        pass
    else:
        tree = Phylo.read(phy_information, 'newick')
        tree_id = select_tree_id(tree)
        tree_id.sort()
        if tree_id != seq_id:
            1/0
except:
    print('Please input correct phy_information parament with the same seq id (newick format)')

#identity_cutoff
try:
    index = sys.argv.index('-identity_cutoff')+1
    identity_cutoff  = sys.argv[index]
except:
    identity_cutoff =0.8

try:
    identity_cutoff = float(identity_cutoff)
    if identity_cutoff > 1 or identity_cutoff <0:
        1/0
except:
    print('Please input correct identity_cutoff (float 0~1, default: 0.8)')
    sys.exit(0)

#density_fold
try:
    index = sys.argv.index('-density_fold')+1
    density_fold  = sys.argv[index]
except:
    density_fold =1.5

try:
    density_fold = float(density_fold)
except:
    print('Please input correct density_fold (float, default: 1.5)')
    sys.exit(0)

#threads
try:
    index = sys.argv.index('-threads')+1
    threads  = sys.argv[index]
except:
    threads = 'AUTO'

try:
    if threads != 'AUTO':
        threads = int(threads)
except:
    print('Please input correct threads (int or AUTO, default: AUTO)')
    sys.exit(0)    

##########################
#clustering new sequences#
##########################
if output_location[-1] != '/':
    output_location = output_location+'/'

if len(combine_location) == 0:
    add_aligned_sequence(add_seq_location = addition_seq_location,
                         aligned_seq_location=aligned_seq_location,
                         output_location=output_location+'combine.fasta',
                         thread=threads)
    combine_location = output_location+'combine.fasta'

new_seqid,out_seq_id,median_list = gap_t_test(add_seq_location = addition_seq_location,
                       aligned_seq_location=aligned_seq_location,
                       IDfold=IDfold,
                       output_location=combine_location)

with open(output_location+'outdatabase.id','w') as f:
    for i in range(len(out_seq_id)):
        f.write(str(out_seq_id[i])+'\n')

pd.DataFrame(median_list).to_csv('./median_list.csv',header=None,index=None)

database = []
with open(database_location +'.match','r') as f:
    for i in f.readlines():
        i = i.split('\n')
        database.append(i[0])
        
if  database[0].split(',')[2] != '':
    convert_rule =  database[0].split(',')[2]
    convert_rule_location = ''
    
    from rpy2.robjects import r as Rcode
    from rpy2.robjects.packages import importr as Rrequire
    Rrequire('ape')
    Rcode("align <- read.dna('%s',format = 'fasta')" % combine_location)

    Rcode("distance <- dist.dna(align,as.matrix = T,model = '%s')" % convert_rule)

    Rcode("write.csv(file = './p_distance_combine.csv',distance)")
    
    pdistance = pd.read_csv('./p_distance_combine.csv')
    
else:
    convert_rule_location = database_location+'.convert'
    pdistance = ''

clustering_result,clustering_result_total,clustering_density,clustering_output = clustering_information_tree(database_location=database_location+'.match',
                                                                                           aligned_seq_location=aligned_seq_location,
                                                                                           combine_seq_location=combine_location,
                                                                                           convert_rule_location=convert_rule_location,
                                                                                           new_seqid=new_seqid,
                                                                                           identity_cutoff=identity_cutoff,
                                                                                           density_fold=density_fold,
                                                                                           pdistance = pdistance)

with open(output_location+'clustering_result.csv','w') as f:
    for i in range(len(clustering_result)):
        f.write(str(clustering_result[i][0])+',')
        f.write(str(clustering_result[i][1])+',')
        f.write(str(clustering_density[i])+'\n')

if  phy_information != '':                                                                                                                                                                                 
    # try:
    combine_tree = add_tree(clustering_result,
                    ML_tree_location=phy_information,
                    aligned_seq_location=aligned_seq_location,
                    combine_seq_location=combine_location,
                    threads=threads
                    )

    Phylo.write(combine_tree,output_location+'combined_tree.nwk','newick')
    # except:
    #     print('Please check the input database and phylogenetic tree')


    
    
# with open('./tree_time_info.csv','w') as f:
#     f.write('total,'+str(total_end-total_start)+'\n')
#     for i in range(len(time_info)):
#         f.write(str(time_info[i][0])+','+str(time_info[i][1])+','+str(time_info[i][2])+'\n')

# np.savetxt('./clustering_output.csv',np.array(clustering_output),delimiter=',')














     
        

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mabin
"""

#####Make Database function
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

def bases_convert (pi, sequence, convert_rule_location = '' ):
    import numpy as np
    from Bio import SeqIO
    if convert_rule_location == '':
        A = np.array([1,0,0,0,1,0])* (1-pi[0])
        C = np.array([0,1,0,0,0,1])* (1-pi[1])
        G = np.array([0,0,1,0,1,0])* (1-pi[2])
        T = np.array([0,0,0,1,0,1])* (1-pi[3])
        
        # A = np.array([1,0,0,0,1,0])
        # C = np.array([0,1,0,0,0,1])
        # G = np.array([0,0,1,0,1,0])
        # T = np.array([0,0,0,1,0,1])
        
        # A = np.array([1,0,0,0])* (1-pi[0])
        # C = np.array([0,1,0,0])* (1-pi[1])
        # G = np.array([0,0,1,0])* (1-pi[2])
        # T = np.array([0,0,0,1])* (1-pi[3])
        
        # A = np.array([1,0,0,0])
        # C = np.array([0,1,0,0])
        # G = np.array([0,0,1,0])
        # T = np.array([0,0,0,1])
        
        # A = np.array([1,0,0,0,1,0])* (pi[0])
        # C = np.array([0,1,0,0,0,1])* (pi[1])
        # G = np.array([0,0,1,0,1,0])* (pi[2])
        # T = np.array([0,0,0,1,0,1])* (pi[3])
        
        # A = np.array([1,0,0,0])* (pi[0])
        # C = np.array([0,1,0,0])* (pi[1])
        # G = np.array([0,0,1,0])* (pi[2])
        # T = np.array([0,0,0,1])* (pi[3])
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
    
    
    seq_change_matrix = []
    
    for i in range(len(sequence)):
        tmp_seq = []
        for j in  range(len(sequence[i])):
            if sequence[i][j] == 'a':
                tmp_seq.append(A)
            
            if sequence[i][j] == 'c':
                tmp_seq.append(C)
            
            if sequence[i][j] == 'g':
                tmp_seq.append(G)
            
            if sequence[i][j] == 't':
                tmp_seq.append(T)
            
            if sequence[i][j] == '-':
                tmp_seq.append(gap)
            
            if sequence[i][j] == 'r':
                tmp_seq.append(R)
                    
            if sequence[i][j] == 'y':
                tmp_seq.append(Y)
        
            if sequence[i][j] == 's':
                tmp_seq.append(S)
    
            if sequence[i][j] == 'w':
                tmp_seq.append(W)
    
            if sequence[i][j] == 'k':
                tmp_seq.append(K)
            
            if sequence[i][j] == 'm':
                tmp_seq.append(M)
    
            if sequence[i][j] == 'b':
                tmp_seq.append(B)            
    
            if sequence[i][j] == 'd':
                tmp_seq.append(D)
    
            if sequence[i][j] == 'h':
                tmp_seq.append(H)
    
            if sequence[i][j] == 'v':
                tmp_seq.append(V)
                
            if sequence[i][j] == 'n':
                tmp_seq.append(N)
        
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
    
def information_clustering(seq_change_matrix_PCA,seq_id,distance_exponent = 2, clustering_method = 'single',clustering_information = '',cluster_number = 2):
    ####make Database
    from sklearn.cluster import AgglomerativeClustering
    from scipy.spatial.distance import pdist, squareform
    import numpy as np 
    import csv
    import pandas as pd
    ####calcuate distance matrix
    if  distance_exponent == 2:
        distance_matrix = pdist(seq_change_matrix_PCA,'euclidean')
        distance_matrix = squareform(distance_matrix)
    elif distance_exponent == 1:
        distance_matrix = pdist(seq_change_matrix_PCA,'cityblock')
        distance_matrix = squareform(distance_matrix)
    else:
        distance_matrix = pdist(seq_change_matrix_PCA,'minkowski',p=distance_exponent)
        distance_matrix = squareform(distance_matrix)
    ####
        
        
    ###clustering
    output_id = []
    output_location = []
    output_identity = []
    output_index = []
    output_density = []
    ### identity = jaccard value
    
    
    if clustering_information == '':
        clustering = AgglomerativeClustering(n_clusters = cluster_number,affinity = 'precomputed',
                                             linkage = clustering_method).fit(distance_matrix)
        for i in range(cluster_number):
            output_id.append('cluster%s' % i)
            output_location.append(np.where(clustering.labels_==i))
            output_identity.append(1)
            output_density.append(np.max(distance_matrix[np.where(clustering.labels_==i)[0],:][:,np.where(clustering.labels_==i)[0]]))
    
    else:
        ###input information 
        information = pd.read_csv(clustering_information, sep=',', header=None)
        ###information -- seq_id, clade, subclade .....
        
        cluster_level_number = len(information.loc[0]) - 1 ##remove seqid
        
        
        seq_id = pd.DataFrame(seq_id)
        information = pd.merge(seq_id,information,on=0) ##match information
        
                
      
        for z in range(1,3):
            if z == 1:
                cluster_information_index = []
                for i in range(len(pd.value_counts(information[z]).index)):
                     #   clustering_number_remove += 1
                    cluster_information_index.append(pd.value_counts(information[z]).index[i])###input information index 
                
                ###identity -> Jaccard A n B/A U B
                
                tmp_cluster_identity = [[] for i in range(len(cluster_information_index))]
                tmp_cluster_location = [[] for i in range(len(cluster_information_index))]
                
                
                
                for clustering_number in range(1,len(cluster_information_index)*3):
                    
                    clustering = AgglomerativeClustering(n_clusters = clustering_number,affinity = 'precomputed',
                                         linkage = clustering_method).fit(distance_matrix)
                
                    for i in range(clustering_number):
                        for j in range(len(pd.value_counts(information[z][list(np.where(clustering.labels_ == i)[0])]).index)):
                            match_information_index = cluster_information_index.index(pd.value_counts(information[z][list(np.where(clustering.labels_ == i)[0])]).index[j])
                            
                            tmp_cluster_map_number = pd.value_counts(information[z][list(np.where(clustering.labels_ == i)[0])])[j]
                            tmp_cluster_total_number = sum(pd.value_counts(information[z][list(np.where(clustering.labels_ == i)[0])]))
                            total_information_number = len(information[information[1] == cluster_information_index[match_information_index]])
                            identity = tmp_cluster_map_number / (tmp_cluster_total_number+total_information_number-tmp_cluster_map_number)
                            
                            tmp_cluster_identity[match_information_index].append(identity)
                            
                            tmp_cluster_location[match_information_index].append(list(np.where(clustering.labels_ == i)[0]))
                            
                    
                for i in range (len(tmp_cluster_identity)):
                    max_identity = max(tmp_cluster_identity[i])
                    max_identity_index = np.where(np.array(tmp_cluster_identity[i]) == max_identity)[0][0]
                    
                    output_id.append(cluster_information_index[i])
                    output_identity.append(max_identity)
                    output_location.append(tmp_cluster_location[i][max_identity_index])
                    output_index.append(z)
                    output_density.append(np.max(distance_matrix[tmp_cluster_location[i][max_identity_index],:][:,tmp_cluster_location[i][max_identity_index]]))
                    
            else:
                clustering_index = z - 1
                for y in range (len(np.where(np.array(output_index)== clustering_index)[0])):
                    ##change distance matrix by output id
                    distance_matrix_change = distance_matrix[output_location[np.where(np.array(output_index)==clustering_index)[0][y]],:][:,output_location[np.where(np.array(output_index)==clustering_index)[0][y]]]
                    
                    information_change = information[z][output_location[np.where(np.array(output_index)==clustering_index)[0][y]]]
                        
                    cluster_information_index = []
                    for i in range(len(pd.value_counts(information_change).index)):
                         #   clustering_number_remove += 1
                        cluster_information_index.append(pd.value_counts(information_change).index[i])###input information index 
                    
                    ###identity -> Jaccard A n B/A U B
                    
                    tmp_cluster_identity = [[] for i in range(len(cluster_information_index))]
                    tmp_cluster_location = [[] for i in range(len(cluster_information_index))]
                    
                    
                    for clustering_number in range(1,len(cluster_information_index)*3):
                        
                        clustering = AgglomerativeClustering(n_clusters = clustering_number,affinity = 'precomputed',
                                             linkage = clustering_method).fit(distance_matrix_change)
                    
                        for i in range(clustering_number):
                            for j in range(len(pd.value_counts(information_change[information_change.index[list(np.where(clustering.labels_ == i)[0])]]).index)):
                                match_information_index = cluster_information_index.index(pd.value_counts(information_change[information_change.index[list(np.where(clustering.labels_ == i)[0])]]).index[j])
                                
                                tmp_cluster_map_number = pd.value_counts(information_change[information_change.index[list(np.where(clustering.labels_ == i)[0])]])[j]
                                tmp_cluster_total_number = sum(pd.value_counts(information_change[information_change.index[list(np.where(clustering.labels_ == i)[0])]]))
                                total_information_number = len(information_change[information_change == cluster_information_index[match_information_index]])
                                identity = tmp_cluster_map_number / (tmp_cluster_total_number+total_information_number-tmp_cluster_map_number)
                                
                                tmp_cluster_identity[match_information_index].append(identity)
                                
                                tmp_cluster_location[match_information_index].append(information_change.index[list(np.where(clustering.labels_ == i)[0])])
                                
                        
                    for i in range (len(tmp_cluster_identity)):
                        max_identity = max(tmp_cluster_identity[i])
                        max_identity_index = np.where(np.array(tmp_cluster_identity[i]) == max_identity)[0][0]
                        
                        output_id.append(cluster_information_index[i])
                        output_identity.append(max_identity)
                        output_location.append(tmp_cluster_location[i][max_identity_index])
                        output_index.append(z)
                        output_density.append(np.max(distance_matrix[tmp_cluster_location[i][max_identity_index],:][:,tmp_cluster_location[i][max_identity_index]]))
    ##################################################
    result = []
    for i in range(len(output_id)):
        result.append([output_id[i],output_location[i],output_index[i],output_identity[i],output_density[i]])
    return result

def ML_tree_clustering(ML_tree_location,seq_change_matrix_PCA,seq_id,max_cluster_number=5,bootstrap_cutoff=90,distance_exponent = 2,clustering_method = 'single'):
    
    from Bio import Phylo
    from Bio.Phylo.BaseTree import Tree
    import pandas as pd
    from sklearn.cluster import AgglomerativeClustering
    from scipy.spatial.distance import pdist, squareform
    import numpy as np 
    ##output
    result = []
    result_index = []
    ##def function 
    
    def select_tree_id(tree):
        tree_id = []
        total_tree_id = tree.get_terminals()
        for i in range(len(total_tree_id)):
            tree_id.append(total_tree_id[i].name)
        return tree_id
        
        ##result output 
        ##ML_location(12122) clustering_location identity density  
        #######
    def calcuate_identity_density(tree_id,seq_id,seq_id_location,distance_matrix,max_cluster_number,clustering_method):
        ###define cluster number
        if distance_matrix.shape[0] < max_cluster_number:
            cluster_number_select = distance_matrix.shape[0]
        else:
            cluster_number_select = max_cluster_number
            
        if seq_id_location == '':
            seq_id_location = list(range(len(seq_id)))
        else:
            seq_id = list(pd.DataFrame(seq_id)[0][seq_id_location])
        ###match cluster result
        result_output = []
        result_max_select = []
        for cluster_number in range(cluster_number_select):
            cluster_number += 1
            clustering = AgglomerativeClustering(n_clusters = cluster_number,affinity = 'precomputed',
             linkage = clustering_method).fit(distance_matrix)
            
            tmp_max_select_list = []
            tmp_result_output = []
            for i in range(cluster_number):
                location_index = np.where(clustering.labels_==i)[0]
                seq_index =  list(pd.DataFrame(seq_id_location)[0][location_index])
                
                clustering_id  = list(pd.DataFrame(seq_id)[0][location_index])
                tmp_inter = len(list(set(tree_id).intersection(set(clustering_id))))
                tmp_union = len(list(set(tree_id).union(set(clustering_id))))
                tmp_identity = tmp_inter / tmp_union
                tmp_density = np.max(distance_matrix[location_index,:][:,location_index])
                
                tmp_max_select_list.append(tmp_identity)
                tmp_result_output.append([location_index,seq_index,tmp_identity,tmp_density])
            
            result_output.append(tmp_result_output[max(np.where(np.array(tmp_max_select_list) == max(tmp_max_select_list))[0])])
            result_max_select.append(max(tmp_max_select_list))
        
        total_result = result_output[max(np.where(np.array(result_max_select) == max(result_max_select))[0])]
       # distance_matrix_change = distance_matrix[total_result[0],:][:,total_result[0]]
        return total_result
    
    def back_search_clade(result,result_index,location,tree_id,seq_id,distance_matrix,clustering_method,max_cluster_number=5):
        location = location[:-1]
        location_list = list(location)
        if len(location_list) == 0 :
            seq_id_location =  ''
            total_result = calcuate_identity_density(tree_id,seq_id,seq_id_location,distance_matrix,max_cluster_number,clustering_method)
       
        elif len(location_list) == 1:
            location_result_index = result_index.index(location)
            seq_id_location = result[location_result_index][1][1]
            distance_matrix1 = distance_matrix[seq_id_location,:][:,seq_id_location]
            total_result1 = calcuate_identity_density(tree_id,seq_id,seq_id_location,distance_matrix1,max_cluster_number,clustering_method)
            
            seq_id_location =  ''
            max_cluster_number = max_cluster_number * 2
            total_result2 = calcuate_identity_density(tree_id,seq_id,seq_id_location,distance_matrix,max_cluster_number,clustering_method)


            if total_result1[2] > total_result2[2]:
                total_result = total_result1
            else:
                total_result = total_result2
                
        else :
            location_result_index1 = result_index.index(location)
            seq_id_location1 = result[location_result_index1][1][1]
            distance_matrix1 = distance_matrix[seq_id_location1,:][:,seq_id_location1]
            total_result1 = calcuate_identity_density(tree_id,seq_id,seq_id_location1,distance_matrix1,max_cluster_number,clustering_method)
            
            location_result_index2 = result_index.index(location[:-1])
            max_cluster_number = max_cluster_number * 2
            seq_id_location2 = result[location_result_index2][1][1]
            distance_matrix2 = distance_matrix[seq_id_location2,:][:,seq_id_location2]
            total_result2 = calcuate_identity_density(tree_id,seq_id,seq_id_location2,distance_matrix2,max_cluster_number,clustering_method)
            
            if total_result1[2] > total_result2[2]:
                total_result = total_result1
            else:
                total_result = total_result2
    
        return total_result
    
    def search_tree(tree,distance_matrix, seq_id,clustering_method,bootstrap_cutoff=90,seq_id_location='',location = ''):
        if len(tree) != 0:
            if str(type(tree.confidence)) != "<class 'NoneType'>":
                if tree.confidence > bootstrap_cutoff:
                    tree_id = select_tree_id(tree)
                    total_result = back_search_clade(result,result_index,location,tree_id,seq_id,distance_matrix,clustering_method,max_cluster_number)
                    #seq_id_location = total_result[1]
                    
                    result.append([location,total_result])
                    result_index.append(location)
                    
                    if len(total_result[0]) != 1:
                        location0 = location + '0'
                        location1 = location + '1'
                        
                        search_tree(tree = tree.clades[0],location = location0,distance_matrix=distance_matrix,clustering_method=clustering_method,
                                    bootstrap_cutoff=90,seq_id=seq_id,seq_id_location=seq_id_location)
                        
                        search_tree(tree = tree.clades[1],location = location1,distance_matrix=distance_matrix,clustering_method=clustering_method,
                                    bootstrap_cutoff=90,seq_id=seq_id,seq_id_location=seq_id_location)
              
    ##input tree
    tree = Phylo.read(ML_tree_location, 'newick')
    ##mid point
    tree.root_at_midpoint()
    
    ####calcuate distance matrix
    if  distance_exponent == 2:
        distance_matrix = pdist(seq_change_matrix_PCA,'euclidean')
        distance_matrix = squareform(distance_matrix)
    elif distance_exponent == 1:
        distance_matrix = pdist(seq_change_matrix_PCA,'cityblock')
        distance_matrix = squareform(distance_matrix)
    else:
        distance_matrix = pdist(seq_change_matrix_PCA,'minkowski',p=distance_exponent)
        distance_matrix = squareform(distance_matrix)
    ####
    search_tree(tree.clade[0],distance_matrix=distance_matrix,seq_id=seq_id,location='0',bootstrap_cutoff=bootstrap_cutoff,clustering_method=clustering_method)
    search_tree(tree.clade[1],distance_matrix=distance_matrix,seq_id=seq_id,location='1',bootstrap_cutoff=bootstrap_cutoff,clustering_method=clustering_method)
    
    return result
        

#####Clustering function
def add_aligned_sequence(add_seq_location, aligned_seq_location,output_location = '',thread = 80):
    import os 
    if output_location == '':
        output_location = aligned_seq_location + '.combine'
        
    os.system("mafft --quiet --thread %s --add %s --reorder %s > %s" %(thread, add_seq_location, aligned_seq_location, output_location))
    
    return 

def gap_t_test(add_seq_location, aligned_seq_location, output_location = '', cutoff_pvalue = 0.01):
    from scipy import stats
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
    
    for i in range(len(add_seq_id)):
        
        add_seq_information = after_information[after_id.index(add_seq_id[i]),:]
        after_test = np.zeros(after_match_information.shape)
        
        for j in range(after_match_information.shape[0]):
            after_test[j,:] =  after_match_information[j,:] - add_seq_information
        
        after_test = sum(after_test.T)
        ##calcuate paired t test
        pvalue = stats.ttest_rel(before_test,after_test)[1]
        
        if pvalue == 'nan':
            new_seqid.append(add_seq_id[i])
            print(add_seq_id[i] + ' IN Database')
        
        else:
            if float(pvalue) >= float(cutoff_pvalue):
                new_seqid.append(add_seq_id[i])
                print(add_seq_id[i] + ' IN Database')
            else:
                print(add_seq_id[i] + ' OUT Database')
    return new_seqid

def clustering_information_tree(database_location,aligned_seq_location,combine_seq_location,new_seqid,convert_rule_location,identity_cutoff=0.9):
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
    distance_exponent = int(database[0].split(',')[1])
    
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
        tmp_search_out_total=[]
        for j in range(len(database)):
            try:
                database[j].split('/')[-1].split(',').index(str(clustering_output[i]))
                tmp_search_out_total.append(database[j].split('/')[0])
                if clustering_density[i] < float(database[j].split('/')[2]):
                    tmp_search_out.append(database[j].split('/')[0])
            except:
                pass
        clustering_DBsearch_id.append(clustering_id[i])
        clustering_DBsearch.append(tmp_search_out)
        clustering_DBsearch_total.append(tmp_search_out_total)
    ######select min-tree 
    clustering_result = []
    for i in range(len(clustering_DBsearch)):
        clustering_result.append([clustering_DBsearch_id[i],clustering_DBsearch[i][-1]])
    
    clustering_result_total = []
    for i in range(len(clustering_DBsearch)):
        clustering_result_total.append([clustering_DBsearch_id[i],clustering_DBsearch_total[i][:]])
    
    return clustering_result

def add_tree(clustering_result,ML_tree_location,aligned_seq_location,combine_seq_location):
    from Bio import SeqIO
    from Bio import Phylo
    from Bio.Phylo.BaseTree import Tree
    import pandas as pd
    import os 
    from Bio.SeqRecord import SeqRecord
    from io import StringIO
    import random
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
    
    for i in reversed(range(len(location_list_index))):
        for j in range(len(location_list_index)):
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
    old_tree = tree
    ##merge same tree location 
    for i in range(len(location_list_index)):
        tree_location =location_list_index[i]
        add_seq_id = add_seq_id_list[i]
        db_seq_id = select_tree_id(read_tree(tree,tree_location))
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
            os.system('iqtree -s ./tmp_out.fasta -T 40 ')
            new_tree = Phylo.read('./tmp_out.fasta.treefile', 'newick')
            new_tree.root_at_midpoint()
            os.system('rm ./tmp_out.*')
        else:
            SeqIO.write(tmp_out_fasta, './tmp_out.fasta', 'fasta')
            os.system('iqtree -s ./tmp_out.fasta -B 1000 -T 40 ')
            new_tree = Phylo.read('./tmp_out.fasta.treefile', 'newick')
            new_tree.root_at_midpoint()
            os.system('rm ./tmp_out.*')
        
        Phylo.write(new_tree.clade,'./tree/%s' % tree_location,'newick')
        
        new_location = './tree/%s' % tree_location+'.csv'
        with open(new_location ,'w') as f:
            for h in range(len(add_seq_id)):
                f.write(str(add_seq_id[h])+'\n')
            
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
        
        print(branch_length)
        
        exec(tree_location_change + '.branch_length= %s' % branch_length)

        
    return tree
        
        
        

    

        
        
    

    
    
        

####################################

import sys
from Bio import Phylo
import numpy as np


#aligned_seq_location = './final_Sgene_align.fasta'
aligned_seq_location = '/Users/mabin/Desktop/phylogenyPCA/cov2021/final/data/final_Sgene_align.fasta'

#ML_tree_location = './final_Sgene_align.fasta.treefile'
ML_tree_location = '/Users/mabin/Desktop/phylogenyPCA/cov2021/final/data/final_Sgene_align.fasta.treefile'

clustering_information='/Users/mabin/Desktop/phylogenyPCA/cov2021/final/data/test_software/addtion_test_information_total_organ.csv'
#addition_seq_location = './S_gene_addition_test.fasta'
addition_seq_location = '/Users/mabin/Desktop/phylogenyPCA/cov2021/final/data/test_software/S_gene_addition_test.fasta'

#combine_seq_location = './final_Sgene_align.combine.fasta'
combine_seq_location = '/Users/mabin/Desktop/phylogenyPCA/cov2021/final/data/test_software/combine.fasta'

database_location = './tree_cluster_result'
#database_location = '/Users/mabin/Desktop/phylogenyPCA/cov2021/final/data/test_software/database'
database_location='/Users/mabin/Downloads/information_cluster_result'

#convert_rule_location ='./convert_matrix.csv'
convert_rule_location = '/Users/mabin/Desktop/phylogenyPCA/cov2021/final/data/test_software/convert_matrix.csv'
#######################
pi,seq_id,sequence= calcuate_bases_frequency(aligned_seq_location)       

seq_change_matrix,convert_matrix= bases_convert(pi,sequence)

seq_change_matrix_PCA = PCA_improved(seq_change_matrix)

result = information_clustering(seq_change_matrix_PCA,seq_id,distance_exponent=2,clustering_method='single',clustering_information=clustering_information)


result = ML_tree_clustering(ML_tree_location,
                            seq_change_matrix_PCA=seq_change_matrix_PCA,seq_id=seq_id)

##write convert matrix
np.savetxt(convert_rule_location,convert_matrix,delimiter=',')
###write tree database
with open(database_location,'w') as f:
    f.write('PCA_on'+'/'+''+','+'2'+'\n')
    for i in range (len(result)):
        tree_location = result[i][0]
        identity = result[i][1][2]
        density = result[i][1][3]
        seq_location = result[i][1][1]
        
        f.write(tree_location+'/'+str(identity)+'/'+str(density)+'/')
        for j in range(len(seq_location)):
            f.write(str(seq_location[j])+',')
        f.write('\n')

#write information database
with open('/Users/mabin/Downloads/information_cluster_result','w') as f:
    f.write('PCA_on'+'/'+''+','+'2'+'\n')
    for i in range (len(result)):
        seq_id = result[i][0]
        identity = result[i][3]
        density = result[i][4]
        seq_location = result[i][1]
        
        f.write(seq_id+'/'+str(identity)+'/'+str(density)+'/')
        for j in range(len(seq_location)):
            f.write(str(seq_location[j])+',')
        f.write('\n')

add_aligned_sequence(addition_seq_location, aligned_seq_location,combine_seq_location,thread=56)


new_seqid = gap_t_test(addition_seq_location,aligned_seq_location,output_location=combine_seq_location,cutoff_pvalue=0)



clustering_result = clustering_information_tree(database_location,aligned_seq_location=aligned_seq_location,combine_seq_location=combine_seq_location,
                                                new_seqid=new_seqid,convert_rule_location=convert_rule_location)
   
tree = add_tree(clustering_result, ML_tree_location, aligned_seq_location, combine_seq_location)

Phylo.write(tree,'./tree_test.nwk','newick')









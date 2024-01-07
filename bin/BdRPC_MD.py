#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
BdRPC_MD.py

                    Bd-RPC 
(Bases dependent Rapid Phylogenetic Clustering)

                MAKE DATABASE
                                
                                Author: Ma Bin
'''
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

def bases_convert(pi, sequence,convert_rule='Method1',convert_rule_location = '' ):
    import numpy as np
    from Bio import SeqIO 
    if convert_rule_location == '' and convert_rule == '':
        A = np.array([1,0,0,0,1,0])* (1-pi[0])
        C = np.array([0,1,0,0,0,1])* (1-pi[1])
        G = np.array([0,0,1,0,1,0])* (1-pi[2])
        T = np.array([0,0,0,1,0,1])* (1-pi[3])
    elif convert_rule_location == '':
        if convert_rule =='Method1':
            A = np.array([1,0,0,0,1,0])* (1-pi[0])
            C = np.array([0,1,0,0,0,1])* (1-pi[1])
            G = np.array([0,0,1,0,1,0])* (1-pi[2])
            T = np.array([0,0,0,1,0,1])* (1-pi[3])
        elif convert_rule=='Method2':
            A = np.array([1,0,0,0])* (1-pi[0])
            C = np.array([0,1,0,0])* (1-pi[1])
            G = np.array([0,0,1,0])* (1-pi[2])
            T = np.array([0,0,0,1])* (1-pi[3])
        elif convert_rule=='Method3':
            A = np.array([1,0,0,0,1,0])* (pi[0])
            C = np.array([0,1,0,0,0,1])* (pi[1])
            G = np.array([0,0,1,0,1,0])* (pi[2])
            T = np.array([0,0,0,1,0,1])* (pi[3])
        elif convert_rule=='Method4':
            A = np.array([1,0,0,0])* (pi[0])
            C = np.array([0,1,0,0])* (pi[1])
            G = np.array([0,0,1,0])* (pi[2])
            T = np.array([0,0,0,1])* (pi[3])
        elif convert_rule=='Method5':
            A = np.array([1,0,0,0,1,0])* 1
            C = np.array([0,1,0,0,0,1])* 1
            G = np.array([0,0,1,0,1,0])* 1
            T = np.array([0,0,0,1,0,1])* 1
        elif convert_rule=='Method6':
            A = np.array([1,0,0,0])* 1
            C = np.array([0,1,0,0])* 1
            G = np.array([0,0,1,0])* 1
            T = np.array([0,0,0,1])* 1
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

def information_clustering(seq_change_matrix_PCA,seq_id,distance_exponent = 2, clustering_method = 'single',clustering_information = '',cluster_number = 2, pdistance = ''):
    ####make Database
    from sklearn.cluster import AgglomerativeClustering
    from scipy.spatial.distance import pdist, squareform
    import numpy as np 
    import pandas as pd
    ####calcuate distance matrix
    if len(pdistance) != 0:
        distance_matrix=np.zeros([len(seq_id),len(seq_id)])
        for i in range(len(seq_id)):
            for j in range(i+1,len(seq_id)):
                index = list(pdistance['Unnamed: 0']).index(seq_id[j])
                distance_matrix[i,j]=pdistance[seq_id[i]][index]
                distance_matrix[j,i]=pdistance[seq_id[i]][index]
    else:
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
            output_index.append(0)
            output_density.append(np.max(distance_matrix[np.where(clustering.labels_==i)[0],:][:,np.where(clustering.labels_==i)[0]]))
    
    else:
        ###input information 
        information = pd.read_csv(clustering_information, sep=',', header=None)
        ###information -- seq_id, clade, subclade .....
        
        cluster_level_number = len(information.loc[0])##remove seqid
        
        
        seq_id = pd.DataFrame(seq_id)
        information = pd.merge(seq_id,information,on=0) ##match information
        
                
        for z in range(1,cluster_level_number):
            if z == 1:
                cluster_information_index = []
                for i in range(len(pd.value_counts(information[z]).index)):
                     #   clustering_number_remove += 1
                    cluster_information_index.append(pd.value_counts(information[z]).index[i])###input information index 
                
                ###Matching Identity -> Jaccard A n B/A U B
                
                tmp_cluster_identity = [[] for i in range(len(cluster_information_index))]
                tmp_cluster_location = [[] for i in range(len(cluster_information_index))]
                
                if len(cluster_information_index)*3 > distance_matrix.shape[0]:
                    max_clustering_number = distance_matrix.shape[0]
                else:
                    max_clustering_number = len(cluster_information_index)*3
                
                for clustering_number in range(1,max_clustering_number):
                    
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
                    
                    if len(cluster_information_index)*3 > distance_matrix_change.shape[0]:
                        max_clustering_number = distance_matrix_change.shape[0]
                    else:
                        max_clustering_number = len(cluster_information_index)*3
                        
                    for clustering_number in range(1,max_clustering_number):
                        
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
                                
                    if tmp_cluster_identity != [[]]:
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

def information_clustering_old(seq_change_matrix_PCA,seq_id,distance_exponent = 2, clustering_method = 'single',clustering_information = '',cluster_number = 2):
    ####make Database
    from sklearn.cluster import AgglomerativeClustering
    from scipy.spatial.distance import pdist, squareform
    import numpy as np 
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
        
        cluster_level_number = len(information.loc[0])##remove seqid
        
        
        seq_id = pd.DataFrame(seq_id)
        information = pd.merge(seq_id,information,on=0) ##match information
        
                
        for z in range(1,cluster_level_number):
            if z == 1:
                cluster_information_index = []
                for i in range(len(pd.value_counts(information[z]).index)):
                     #   clustering_number_remove += 1
                    cluster_information_index.append(pd.value_counts(information[z]).index[i])###input information index 
                
                ###Matching Identity -> Jaccard A n B/A U B
                
                tmp_cluster_identity = [[] for i in range(len(cluster_information_index))]
                tmp_cluster_location = [[] for i in range(len(cluster_information_index))]
                
                if len(cluster_information_index)*3 > distance_matrix.shape[0]:
                    max_clustering_number = distance_matrix.shape[0]
                else:
                    max_clustering_number = len(cluster_information_index)*3
                
                for clustering_number in range(1,max_clustering_number):
                    
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
                    
                    if len(cluster_information_index)*3 > distance_matrix_change.shape[0]:
                        max_clustering_number = distance_matrix_change.shape[0]
                    else:
                        max_clustering_number = len(cluster_information_index)*3
                        
                    for clustering_number in range(1,max_clustering_number):
                        
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
                                
                    if tmp_cluster_identity != [[]]:
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

def select_tree_id(tree):
    tree_id = []
    total_tree_id = tree.get_terminals()
    for i in range(len(total_tree_id)):
        tree_id.append(total_tree_id[i].name)
    return tree_id
        
def ML_tree_clustering(ML_tree_location,seq_change_matrix_PCA,seq_id,max_cluster_number=5,bootstrap_cutoff=90,distance_exponent = 2,clustering_method = 'single',pdistance = ''):
    
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


            if float(total_result1[2]) > float(total_result2[2]):
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
            
            if float(total_result1[2]) > float(total_result2[2]):
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
                                    bootstrap_cutoff=bootstrap_cutoff,seq_id=seq_id,seq_id_location=seq_id_location)
                        
                        search_tree(tree = tree.clades[1],location = location1,distance_matrix=distance_matrix,clustering_method=clustering_method,
                                    bootstrap_cutoff=bootstrap_cutoff,seq_id=seq_id,seq_id_location=seq_id_location)
              
    ##input tree
    tree = Phylo.read(ML_tree_location, 'newick')
    ##mid point
    tree.root_at_midpoint()
    
    ####calcuate distance matrix    
    if len(pdistance) != 0:
        distance_matrix=np.zeros([len(seq_id),len(seq_id)])
        for i in range(len(seq_id)):
            for j in range(i+1,len(seq_id)):
                index = list(pdistance['Unnamed: 0']).index(seq_id[j])
                distance_matrix[i,j]=pdistance[seq_id[i]][index]
                distance_matrix[j,i]=pdistance[seq_id[i]][index]
    else:
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
    

    for i in range(len(result)):
        result[i] = [result[i][0],result[i][1][1],0,result[i][1][2],result[i][1][3]]
    return result
        

###############################

import sys
import os
from Bio import SeqIO
from Bio import Phylo
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-align',help="Location of aligned sequences. (required) [no punctuation mark: '/' or ',']")
parser.add_argument('-o',help="Directory to store the result. (required)")
parser.add_argument('-convert_method',help="Select the recoding method from Method 1 to Method 6, and the distance estimated methods designed by the R package('ape') [raw, N, TS, TV, JC69, K80,F81, K81, F84, BH87, T92, TN93, GG95, logdet, paralin, indel, indelblock] (default: 'Method1')")
parser.add_argument('-seq_convert',help="Location of convert matrix, the script will use (1-pi,0,0,0,1-pi,0) as default (Method 1).")
parser.add_argument('-PCA',help="Use PCA program to increase the speed or not [on or off]. (default: 'on')")
parser.add_argument('-PCAcomponents',help="If '-PCA' is on, '-PCAcomponents' can be set as the PCA components. (<=number of the sequences and <= length of recoding sequences) (default: max)")
parser.add_argument('-dis_exponent',help="The exponent of minkowski distance. (default: 2)")
parser.add_argument('-Cmethod',help="The method of hierarchical clustering. (single, average, complete, ward) (default: single)")
parser.add_argument('-tax_information',help="The location of sequences taxonomy information. (csv file) [seq_id,clade,subclade,sub-subclade....] [no punctuation mark: '/' or ',']")
parser.add_argument('-phy_information',help="The location of tree with newick format. [no punctuation mark: '/' or ',']")
parser.add_argument('-Cnumber',help="If '-tax_information and -phy_information' not apply, the numebr of cluster will be calcuated without identity. (default: 5)")
parser.add_argument('-bootstrap_cutoff',help="The cutoff value to stop the tree traversal. (default: 90)")
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
    print('Please input aligned sequences')
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

#convert method
try:
    index = sys.argv.index('-convert_method')+1
    convert_rule  = sys.argv[index]
except:
    convert_rule = ''

if convert_rule == '':
    pass
else:
    if convert_rule not in ['Method1','Method2','Method3','Method4','Method5','Method6',
                            "raw", "N", "TS", "TV", "JC69", "K80","F81", "K81", "F84", 
                            "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel",
                            "indelblock"]:
        print('Please input correct recoing method')
        sys.exit(0)

#convert matrix
try:
    index = sys.argv.index('-seq_convert')+1
    convert_rule_location  = sys.argv[index]
except:
    convert_rule_location = ''
    
try:
    if convert_rule_location == '':
        pass
    else:
        # convert_rule = np.loadtxt(convert_rule_location ,delimiter = ',',encoding = 'utf-8-sig') ###sort by A C G T
        if convert_rule.shape[0] != 4:
            1/0
except:
    print('Please input correct convert matrix (CSV format 4*N)')
    sys.exit(0)

#PCA
try:
    index = sys.argv.index('-PCA')+1
    PCA  = sys.argv[index]
except:
    PCA = 'on'
    
try:
    if PCA not in ['on','off']:
        1/0
except:
    print('Please input correct PCA parament (on / off)')
    sys.exit(0)

#PCAcomponents
try:
    index = sys.argv.index('-PCAcomponents')+1
    PCAcomponents  = sys.argv[index]
except:
    PCAcomponents = 'max'

if convert_rule_location=='':
    if seq_number > 6*seq_length:
        PCA = 'off'
else:
    if seq_number > convert_rule.shape[1]*seq_length:
        PCA = 'off'
        
try:
    if PCAcomponents == 'max':
        pass
    elif PCA == 'off':
        pass
    else:
        if int(PCAcomponents) < seq_number:
            1/0
except:
    print('Please input correct PCAcomponents parament (>= number of the sequences)')
    sys.exit(0)

#dis_exponent
try:
    index = sys.argv.index('-dis_exponent')+1
    dis_exponent  = sys.argv[index]
except:
    dis_exponent = 2

try:
    dis_exponent = float(dis_exponent)
except:
    print('Please input correct dis_exponent parament (float)')
    sys.exit(0)

#clustering_method
try:
    index = sys.argv.index('-Cmethod')+1
    clustering_method  = sys.argv[index]
except:
    clustering_method = 'single'

try:
    if clustering_method not in ['single','average','complete','ward']:
        1/0
except:
    print('Please input correct clustering_method parament (single, average, complete, ward)')
    sys.exit(0)

#tax_information
try:
    index = sys.argv.index('-tax_information')+1
    tax_information  = sys.argv[index]
except:
    tax_information = ''


try:
    if tax_information != '':
        tax_sort = list(pd.read_csv(tax_information,header=None)[0])
        tax_sort.sort()
        if tax_sort != seq_id:
            print(tax_sort)
            1/0
        else:
            pd.read_csv(tax_information,header=None)[1]
except:
    print('Please input correct tax_information parament with the same seq id(csv file) [seq_id,clade,subclade,sub-subclade....]')
    sys.exit(0)       
#phy_information
try:
    index = sys.argv.index('-phy_information')+1
    phy_information  = sys.argv[index]
except:
    phy_information = ''

try:
    if tax_information != '':
        print('The tax_information is not empty. The tax_information will be chosen for database creation.')
    elif phy_information == '':
        pass
    else:
        tree = Phylo.read(phy_information, 'newick')
        tree_id = select_tree_id(tree)
        tree_id.sort()
        if tree_id != seq_id:
            1/0
except:
    print('Please input correct phy_information parament with the same seq id (newick format)')
    sys.exit(0)


#clustering_number
try:
    index = sys.argv.index('-Cnumber')+1
    clustering_number  = sys.argv[index]
except:
    clustering_number = 5

try:
    if tax_information != '':
        print('The tax_information is not empty. The tax_information will be chosen for database creation.')
    elif phy_information != '':
        print('The phy_information is not empty. The phy_information will be chosen for database creation.')
    else:
        int(clustering_number)
except:
    print('Please input correct clustering_number parament (int)')
    sys.exit(0)
    
#bootstrap_cutoff
try:
    index = sys.argv.index('-bootstrap_cutoff')+1
    bootstrap_cutoff  = sys.argv[index]
except:
    bootstrap_cutoff = 90


try:
    bootstrap_cutoff = float(bootstrap_cutoff)
    if bootstrap_cutoff > 100 or bootstrap_cutoff < 0:
        1/0
except:
    print('Please input correct bootstrap_cutoff parament (0~100)')
    sys.exit(0)

####################
#establish database#
####################

##tmp test
# aligned_seq_location= '/Users/mabin/Documents/bd-rpc/example/S_gene_align.fasta'
# convert_rule_location = ''
# PCAcomponents='max'
# dis_exponent = 2
# clustering_method='single'
# phy_information = '/Users/mabin/Documents/bd-rpc/example/S_gene.nwk'
# phy_information = ''
# clustering_number=10
# bootstrap_cutoff=90
# tax_information = '/Users/mabin/Documents/bd-rpc/example/S_gene_taxonomy.csv'
# PCA='on'
# convert_rule = 'K80'

pi,seq_id,sequence= calcuate_bases_frequency(aligned_seq_location)       

if convert_rule == '':  
    seq_change_matrix,convert_matrix= bases_convert(pi,sequence,convert_rule_location= convert_rule_location,convert_rule=convert_rule)
    if PCA == 'on':
        seq_change_matrix_PCA = PCA_improved(seq_change_matrix,PCA_components=PCAcomponents)
    else:
        seq_change_matrix_PCA = seq_change_matrix
    pdistance=''
    convert_rule=''
elif len(convert_rule.split('thod')) == 2:
    seq_change_matrix,convert_matrix= bases_convert(pi,sequence,convert_rule_location= convert_rule_location,convert_rule=convert_rule)
    if PCA == 'on':
        seq_change_matrix_PCA = PCA_improved(seq_change_matrix,PCA_components=PCAcomponents)
    else:
        seq_change_matrix_PCA = seq_change_matrix
    pdistance=''
    convert_rule= ''
else:
    from rpy2.robjects import r as Rcode
    from rpy2.robjects.packages import importr as Rrequire
    Rrequire('ape')

    Rcode("align <- read.dna('%s',format = 'fasta')" % aligned_seq_location)

    Rcode("distance <- dist.dna(align,as.matrix = T,model = '%s')" % convert_rule)

    Rcode("write.csv(file = './p_distance.csv',distance)")
    
    pdistance = pd.read_csv('./p_distance.csv')
    
    seq_change_matrix_PCA = ''


if tax_information != '':    
    result = information_clustering(seq_change_matrix_PCA,
                                    seq_id,distance_exponent=dis_exponent,
                                    clustering_method=clustering_method,
                                    clustering_information=tax_information,
                                    cluster_number = clustering_number,
                                    pdistance=pdistance)
else:
    if phy_information=='':
        result = information_clustering(seq_change_matrix_PCA,
                                    seq_id,distance_exponent=dis_exponent,
                                    clustering_method=clustering_method,
                                    clustering_information=tax_information,
                                    cluster_number = clustering_number,
                                    pdistance=pdistance)
    else:
        result = ML_tree_clustering(ML_tree_location=phy_information,
                                    seq_change_matrix_PCA=seq_change_matrix_PCA,
                                    seq_id=seq_id,
                                    bootstrap_cutoff=bootstrap_cutoff,
                                    distance_exponent=dis_exponent,
                                    pdistance = pdistance)

##output database
if output_location[-1] != '/':
    output_location = output_location+'/'
with open(output_location+'database.match','w') as f:
    f.write(str(PCA)+'/'+str(PCAcomponents)+','+str(dis_exponent)+','+convert_rule+'\n')
    for i in range (len(result)):
        seq_id = result[i][0]
        identity = result[i][3]
        density = result[i][4]
        seq_location = result[i][1]
        
        f.write(seq_id+'/'+str(identity)+'/'+str(density)+'/')
        for j in range(len(seq_location)):
            f.write(str(seq_location[j])+',')
        f.write('\n')

if convert_rule != '' and len(convert_rule.split('thod')) == 1:
    pass
else:
    np.savetxt(output_location+'database.convert',convert_matrix,delimiter=',')

        











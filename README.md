# Bases dependent Rapid Phylogenetic Clustering <br>(Bd-RPC)
* Bases dependent Rapid Phylogenetic Clustering (Bd-RPC) is an efficient tool for phylogenetic cluster identification and new sample placement. The Bd-RPC contains two major modules (Make Database and Clustering New Sequences). <br><br>
* In the Make Database module, the aligned sequences will be recoded following the convert matrix, and the components of recoded aligned sequences’ matrix can be extracted for accelerating distance calculation with the principal component analysis (PCA). Then, the distance among sequences will be estimated by the recoded sequence matrix using the Minkowski Distance and the distance matrix will be chosen to match the background information (taxonomy information or phylogenetic tree) for database creation by various Hierarchical Clustering methods and the Simulated Annealing Search algorithm.<br><br>
* After the Bd-RPC database establishment, the new sequences will be added into the aligned sequences using MAFFT and the indel characters will be counted for foreign sequences’ recognition to distinguish whether the new sequences belong to the database. Then, the remainder sequences will be classified according to the Bd-RPC database through the Matching Identity cutoff and clusters’ density. Finally, for the phylogenetic database, the new sequences can be placed into the phylogenetic tree based on the clustering results.<br><br>
* The online toolkit is available on [www.bd-rpc.xyz](http://www.bd-rpc.xyz)

## Installation
### OS Requirements
This package is supported for macOS and Linux. The package has been tested on the following systems:
* macOS: Mojave (10.14.1)
* Linux: Ubuntu (18.04.5)
### Dependencies
Python 3+ 
* numpy
* scipy
* pandas
* biopython
* scikit-learn
* rpy2


R
* ape

Base
* IQ-TREE2 (necessary in phylogenetic tree construction)
* MAFFT

<br>

If you're having difficulties constructing the essential scientific Python packages, we recommend using the [conda](https://docs.conda.io/en/latest/miniconda.html) package/environment manager. <br>

    conda create -n bd_rpc python=3
    conda activate bd_rpc
    conda install numpy scipy pandas biopython scikit-learn rpy2
    conda install -c bioconda mafft iqtree
### Download
    git clone https://github.com/Bin-Ma/bd-rpc.git
    cd bd-rpc/bin

## Manual
This program can recode the aligned sequences to a list of number and match to the background information or phylogenetic tree through hierarchical clustering. For increasing the speed of this program, PCA improvement module can be selected for calcuating the distance between sequences.<br>
### Part 1 -- Make Database
BdRPC_MD.py<br><br>
Usage:

    BdRPC_MD.py [options] -align <location> -o <location>
Output: Bd-RPC database [cluster_location/identity/density/seq_location]
<br>

| Basic options:    |                                                                                                                                                    | 
|:------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------|
| -h                | help                                                                                                                                               |
| -align            | Location of aligned sequences. (required) [no punctuation mark: '/' or ',']                                                                        |
| -o                | Directory to store the result. (required)                                                                                                          |
| -convert_method   | Distance calculation methods (six recoded distance metrics, and seventeen uncoded genetic distance metrics from R package "ape")                   |
| -seq_convert      | Location of convert matrix, the script will use (1-pi,0,0,0,1-pi,0) as default (method 1).                                                         |
| -PCA [on or off]  | Use PCA program to increase the speed or not. (default: 'on')                                                                                      |
| -PCAcomponents    | If "-PCA" is on, '-PCAcomponents' can be set as the PCA components. (<=number of the sequences and <= length of recoding sequences) (default: max) |
| -dis_exponent     | The exponent of minkowski distance. (default: 2)                                                                                                   |
| -Cmethod          | The method of hierarchical clustering. (single, average, complete, ward) (default: single)                                                         |
| -tax_information  | The location of sequences taxonomy information. (csv file) [seq_id,clade,subclade,sub-subclade....] [no punctuation mark: '/' or ',']              |
| -phy_information  | The location of tree with newick format. [no punctuation mark: '/' or ',']                                                                         |
| -Cnumber          | If '-tax_information and -phy_information' not apply, the numebr of cluster will be calcuated without identity. (default: 5)                       |
| -bootstrap_cutoff | The cutoff value to stop the tree traversal. (default: 90)                                                                                         |
<br>

### Part 2 -- Clustering new sequences
BdRPC_CNS.py<br><br>
Usage:

    BdRPC_CNS.py [options] -align <location> -new <location> -o <location> -db <location>

Output: gap-t-test result [seq_id in/out] / clustering result [seq_id cluster_name/tree_location] / combined tree
<br>

| Basic options:   |                                                                                                                          |
|:-----------------|:-------------------------------------------------------------------------------------------------------------------------|
| -align           | Location of aligned sequences. (required) [no punctuation mark: '/' or ',']                                              |
| -new             | Location of new sequences. (required)                                                                                    |
| -o               | Directory to store the result. (required)                                                                                |
| -db              | Location of Bd-RPC database. (required)                                                                                  |
| -IDfold          | The fold of median value in Indel Test. (default: 1.1)                                                                   |
| -phy_information | Location of phylogentics tree. (if the tree is available, the new sequences will be inserted into the phylogenetic tree) |
| -identity_cutoff | The cutoff value of clusters' identity. (0~1, default: 0.8)                                                              |
| -density_fold    | The fold of clusters' density for new samples clustering. (default: 1.5)                                                 |
| -threads         | Threads of mafft align and iqtree. (int, default: auto)                                                                  |

### Part 3 -- BdRPCpackage
This is a python package used to perform phylogenetic new sample placement. Users can flexibly combine functions to achieve personalized functionality.
##### Installation
    pip install BdRPCpackage
##### Function
    BdRPCpackage.calcuate_bases_frequency(aligned_seq_location)
| Input                 | Description                | Output                                                   |
|:----------------------|:---------------------------|:---------------------------------------------------------|
| aligned_seq_location  | Aligned sequence location. | A/C/G/T frequency; Aligned sequence id; Aligned Sequence |

<br>

    BdRPCpackage.bases_convert(pi, sequence,convert_rule='Method1',convert_rule_location = '' )

| Input                  | Description                                  | Output                              |
|:-----------------------|:---------------------------------------------|:------------------------------------|
| pi                     | A/C/G/T frequency.                           | Convert result of aligned sequences |
| sequence               | Aligned sequence.                            | Recoding Conver Rule matrix         |
| convert_rule           | Recoding Convert Rule (select Method 1~6).   |                                     |
| convert_rule_location  | Recoding Convert Rule (A/C/G/T, csv format). |                                     |
<br>

    BdRPCpackage.PCA_improved(seq_change_matrix,PCA_components = 'max')
    #If you do not want to perform PCA on the convert result, skip this function.

| Input                     | Description                          | Output      |
|:--------------------------|:-------------------------------------|:------------|
| seq_change_matrix         | Convert result of aligned sequences. | PCA result  |
| PCA_components            | PCA components number.               |             |
<br>

    BdRPCpackage.information_clustering(seq_change_matrix_PCA,seq_id,distance_exponent = 2, clustering_method = 'single',clustering_information = '',cluster_number = 2, pdistance = '')

| Input                      | Description                                                                                               | Output                                    |
|:---------------------------|:----------------------------------------------------------------------------------------------------------|:------------------------------------------|
| seq_change_matrix_PCA      | Convert result of aligned sequences (perform PCA or not).                                                 | Clustering Result containing:             |
| seq_id                     | Aligned sequence id.                                                                                      | Cluster name (taxonomy)                   |
| distance_exponent          | Distance exponent used in distance calculation among sequences.                                           | Sequences index of Cluster                |
| clustering_method          | Hierarchical clustering methods used to match taxonomy information and sequence recoding distance matrix. | Cluster index (equal to the column index) |
| clustering_information     | Location of taxonomy information (csv format).                                                            | Cluster identity                          |
| cluster_number             | Clustering number of sequence recoding distance matrix (If taxonomy information is empty).                | Cluster density                           |
| pdistance                  | Other distance matrix used to instead sequence recoding distance matrix.                                  |                                           |
<br>

    BdRPCpackage.ML_tree_clustering(ML_tree_location,seq_change_matrix_PCA,seq_id,max_cluster_number=5,bootstrap_cutoff=90,distance_exponent = 2,clustering_method = 'single',pdistance = '')

| Input                    | Description                                                                                                        | Output                                 |
|:-------------------------|:-------------------------------------------------------------------------------------------------------------------|:---------------------------------------|
| ML_tree_location         | Phylogenetic tree location.                                                                                        | Clustering Result containing:          |
| seq_change_matrix_PCA    | Convert result of aligned sequences (perform PCA or not).                                                          | Cluster name (phylogenetic tree node), |
| seq_id                   | Aligned sequence id.                                                                                               | Sequences index of Cluster,            |
| distance_exponent        | Distance exponent used in distance calculation among sequences.                                                    | Cluster index (equal to 0),            |
| clustering_method        | Hierarchical clustering methods used to match taxonomy information and sequence recoding distance matrix.          | Cluster identity,                      |
| max_cluster_number       | Max Hierarchical clustering number in each node ( >=bootstrap cutoff).                                             | Cluster density.                       |
| bootstrap_cutoff         | Bootstrap cutoff in tree searching. If the bootstrap is less than the bootstrap cutoff, the tree search will stop. |                                        |
| pdistance                | Other distance matrix used to instead sequence recoding distance matrix.                                           |                                        |
<br>

    BdRPCpackage.add_aligned_sequence(add_seq_location, aligned_seq_location,output_location = '',thread = 1)
    #This function will combine new sequence and aligned sequences to a multiple sequence alignment using MAFFT.
| Input                      | Description                 | Output                                |
|:---------------------------|:----------------------------|:--------------------------------------|
| add_seq_location           | New sequence location.      | Combined multiple sequence alignment. |
| aligned_seq_location       | Aligned sequence location.  |                                       |
| output_location            | Combined sequence location. |                                       |
| thread                     | Threads of the MAFFT        |                                       |
<br>

    BdRPCpackage.gap_t_test(add_seq_location, aligned_seq_location, output_location = '',IDfold=1.1)
    #This function will combine new sequence and aligned sequences to a multiple sequence alignment using MAFFT.
| Input                      | Description                                              | Output                                                  |
|:---------------------------|:---------------------------------------------------------|:--------------------------------------------------------|
| add_seq_location           | New sequence location.                                   | New sequence ID belongs to the aligned sequence.        |
| aligned_seq_location       | Aligned sequence location.                               | New sequence ID doesn't belong to the aligned sequence. |
| output_location            | Combined sequence location.                              |                                                         |
| IDfold                     | Indel fold change (used to recognize foreign sequences). |                                                         |
<br>

    BdRPCpackage.clustering_information_tree(database_location,aligned_seq_location,combine_seq_location,new_seqid,convert_rule_location,identity_cutoff=0.8,density_fold=1.5,pdistance='',convert_rule='Method1')
| Input                      | Description                                                               | Output                                     |
|:---------------------------|:--------------------------------------------------------------------------|:-------------------------------------------|
| database_location          | BdRPC database location                                                   | New Sequence Clustering Result containing: |
| aligned_seq_location       | Aligned sequence location.                                                | New sequence id,                           |
| combine_seq_location       | Combined sequence location.                                               | New sequnece's cluster.                    |
| new_seqid                  | New sequence ID belongs to the aligned sequence.                          |                                            |
| convert_rule_location      | Location of Recoding Conver Rule matrix.                                  | Min distance to the aligned sequences.     |
| identity_cutoff            | Identity cutoff of cluster searching.                                     |                                            |
| density_fold               | Max fold change of each cluster density.                                  |                                            |
| pdistance                  | Other distance matrix used to instead sequence recoding distance matrix.  |                                            |
| convert_rule               | Recoding Convert Rule (select Method 1~6).                                |                                            |
<br>

    BdRPCpackage.add_tree(clustering_result,ML_tree_location,aligned_seq_location,combine_seq_location,threads=1)
| Input                 | Description                    | Output                      |
|:----------------------|:-------------------------------|:----------------------------|
| clustering_result     | New Sequence Clustering Result | Combined phylogenetic tree. |
| ML_tree_location      | Phylogenetic tree location.    |                             |
| combine_seq_location  | Combined sequence location.    |                             |
| threads               | Threads of the IQ-TREE2        |                             |
<br>

## Example
### _Betacoronavirus_
The example data is located in './example/' folder.
<br>

#### Make Database 
    cd bd-rpc

    #Taxonomy database
    ./bin/BdRPC_MD.py -align ./example/S_gene_align.fasta -o ./example/ -tax_information ./example/S_gene_taxonomy.csv

    #Phylogentic database
    ./bin/BdRPC_MD.py -align ./example/S_gene_align.fasta -o ./example/ -phy_information ./example/S_gene.nwk
<br>

#### Clustering New Sequence
    #Without phylogenetic tree
    ./bin/BdRPC_CNS.py -align ./example/S_gene_align.fasta -new ./example/S_gene_new.fasta -db ./example/database -o ./example/

    #With phylogenetic tree
    ./bin/BdRPC_CNS.py -align ./example/S_gene_align.fasta -new ./example/S_gene_new.fasta -db ./example/database -o ./example/ -phy_information ./example/S_gene.nwk
<br>

### _Alphacoronavirus_
The example data is located in './example/Alpha_Cov.zip'.
<br>

#### Make Database 
    cd bd-rpc
    unzip ./example/Alpha_CoV.zip

    #Phylogenetic database
    ./bin/BdRPC_MD.py ./BdRPC_MD.py -align ./example/ORF1ab_ref.fasta -phy_information ./example/ORF1ab_high_align.fasta.treefile   -o ./
<br>

#### Clustering New Sequence
    #With phylogenetic tree
    ./bin/BdRPC_CNS.py -align ./example/ORF1ab_ref.fasta -new ./example/ORF1ab_test.fasta  -o ./ -db ./database -IDfold 2  -phy_information ./example/ORF1ab_high_align.fasta.treefile
<br>

###  _Alphaherpesvirinae_
The example data is located in './example/Alpha_Cov.zip'.
<br>

#### Make Database 
    cd bd-rpc
    unzip ./example/Alpha_her.zip

    #Phylogenetic database
    ./bin/BdRPC_MD.py ./BdRPC_MD.py -align ./example/her_ref.fasta -phy_information ./example/alpha_her_high_root.treefile   -o ./
<br>

#### Clustering New Sequence
    #With phylogenetic tree
    ./bin/BdRPC_CNS.py -align ./example/her_ref.fasta -new ./example/her_test.fasta  -o ./ -db ./database -IDfold 2  -phy_information ./example/alpha_her_high_root.treefile -convert_method logdet
<br>

### Results
* **database.match** <br>
    Database information containing creation information and matching result.<br><br>
* **database.convert** <br>
    Convert rule using in database creation.<br><br>
* **combine.fasta** <br>
    Combined multiple sequence alignment generated using MAFFT.<br><br>
* **outdatabase.id** <br>
    Sequence ID recognized as the foreign sequence using Indel Recognition.<br><br>
* **clustering_result.csv** <br>
    Column 1 : new sample ID<br>
    Column 2 : cluster name <br>
    Column 3 : min distance to the cluster<br><br>
* **combined_tree.nwk** <br>
    Combined phylogenetic tree with new samples.

## Window system
Users need to set the R into the environment variables. If the R package "ape" is available, users can select 17 uncoded genetic distance. Here, MAFFT program is necessary for the new sample placement and the IQ-TREE2 is used for phylogenetic tree construction, which users should set into the environment variables previously.
### linux
    ./Bd-RPC
### Mac
    Double click or ./Bd-RPC
### Windows
    Double click or open in Rstudio




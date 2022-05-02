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
### Python dependencies
Python 3+<br>

* numpy
* scipy
* pandas
* biopython
* scikit-learn
* csv
<br>

If you're having difficulties constructing the essential scientific Python packages, we recommend using the [conda](https://docs.conda.io/en/latest/miniconda.html) package/environment manager. <br>

    conda create -n bd_rpc python=3
    conda activate bd_rpc
    conda install numpy scipy pandas biopython scikit-learn csv
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

| Basic options: |   | 
| :-----| :---- |
| -align | Location of aligned sequences. (required) [no punctuation mark: '/' or ','] |
| -o | Directory to store the result. (required) |
| -seq_convert |Location of convert matrix, the script will use (1-pi,0,0,0,1-pi,0) as default (method 1).|
|-PCA [on or off]| Use PCA program to increase the speed or not. (default: 'on')|
|-PCAcomponents|If "-PCA" is on, '-PCAcomponents' can be set as the PCA components. (<=number of the sequences and <= length of recoding sequences) (default: max)|
|-dis_exponent|The exponent of minkowski distance. (default: 2)|
|-Cmethod|The method of hierarchical clustering. (single, average, complete, ward) (default: single)|
|-tax_information|The location of sequences taxonomy information. (csv file) [seq_id,clade,subclade,sub-subclade....] [no punctuation mark: '/' or ',']|
|-phy_information|The location of tree with newick format. [no punctuation mark: '/' or ',']|
|-Cnumber|If '-tax_information and -phy_information' not apply, the numebr of cluster will be calcuated without identity. (default: 5)|
|-bootstrap_cutoff|The cutoff value to stop the tree traversal. (default: 90)|
<br>
### Part 2 -- Clustering new sequences
BdRPC_CNS.py<br><br>
Usage:

    BdRPC_CNS.py [options] -align <location> -new <location> -o <location> -db <location>

Output: gap-t-test result [seq_id in/out] / clustering result [seq_id cluster_name/tree_location] / combined tree
<br>

| Basic options: |  |
| :-----| :---- |
| -align  | Location of aligned sequences. (required) [no punctuation mark: '/' or ','] |
| -new  | Location of new sequences. (required) |
|-o|Directory to store the result. (required)|
|-db|Location of Bd-RPC database. (required)|
|-IDfold|The fold of median value in Indel Test. (default: 1.1)|
|-phy_information|Location of phylogentics tree. (if the tree is available, the new sequences will be inserted into the phylogenetic tree)|
|-identity_cutoff|The cutoff value of clusters' identity. (0~1, default: 0.8)|
|-density_fold|The fold of clusters' density for new samples clustering. (default: 1.5)|
|-threads|Threads of mafft align and iqtree. (int, default: 1)|








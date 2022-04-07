# Bases dependent Rapid Phylogenetic Clustering (Bd-RPC)
* Bases dependent Rapid Phylogenetic Clustering (Bd-RPC) is an efficient tool for phylogenetic cluster identification and new sample placement. The Bd-RPC contains two major modules (Make Database and Clustering New Sequences). <br><br>
* In the Make Database module, the aligned sequences will be recoded following the convert matrix, and the components of recoded aligned sequences’ matrix can be extracted for accelerating distance calculation with the principal component analysis (PCA). Then, the distance among sequences will be estimated by the recoded sequence matrix using the Minkowski Distance and the distance matrix will be chosen to match the background information (taxonomy information or phylogenetic tree) for database creation by various Hierarchical Clustering methods and the Simulated Annealing Search algorithm.<br><br>
* After the Bd-RPC database establishment, the new sequences will be added into the aligned sequences using MAFFT and the indel characters will be counted for foreign sequences’ recognition to distinguish whether the new sequences belong to the database. Then, the remainder sequences will be classified according to the Bd-RPC database through the Matching Identity cutoff and clusters’ density. Finally, for the phylogenetic database, the new sequences can be placed into the phylogenetic tree based on the clustering results.<br>

## Installation






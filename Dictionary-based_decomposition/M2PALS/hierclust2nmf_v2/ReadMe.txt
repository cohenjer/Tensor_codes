hierclust2nmf.m is the main function, performing the hierarchical clustering of hyperspectral images. (Note that is can also be applied to any other kind of data.) 

Run the file RunMe.m to have an example with the Urban dataset. 



Modifications w.r.t. v.1: 
- Fix a bug when the input matrix had repeated columns. 
- New output: the index set of the columns of the input matrix correpsonding to the endmembers (=cluster centroids). 
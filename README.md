The script used for analyses in "Novel metric for hyperbolic phylogenetic tree embeddings"

## Reference

## Requirements

The following libraries are used.

* hydra
* MASS
* ape
* ROCR
* RSpectra

## Full node embeddings

```
Rscript embedding_Dall.R
```

The MSE values of each embedding method for the i-th simulation tree are saved in "out/MSE_Dall_1.csv".
The mean MSE values are saved in "out/MSE_Dall_all.csv".

## External node embeddings and folding-in internal nodes

```
Rscript embedding_Dall.R
```

The MSE values of each embedding method for the i-th simulation tree are saved in "out/MSE_Dint_i.csv".
The mean MSE values are saved in "out/MSE_Dint_all.csv".

## External node embeddings and folding-in internal nodes of the partial tree

```
Rscript embedding_Dall.R
```

The MSE values of each embedding method for the i-th simulation tree are saved in "out/MSE_Dpart_i.csv".
The mean MSE values are saved in "out/MSE_Dpart_all.csv".

## Prediction of the MRCA

```
Rscript predicting_MRCA.R
```

The AUC values of each embedding method for the i-th simulation tree are saved in "out/AUC_MRCA_i.csv".
The mean AUC values are saved in "out/AUC_MRCA_all.csv".

## Prediction of nearst-neighbor nodes in the partial tree

```
Rscript predicting_nearest_neighbor.R
```

The rank of each embedding method for the i-th simulation tree are saved in "out/Result_NN_MDS_i.csv ", "out/Result_NN_H1_i.csv", and "out/Result_NN_H2_i.csv", respectively.
The mean rank are saved in "out/Result_NN_MDS_all.csv ", "out/Result_NN_H1_all.csv", and "out/Result_NN_H2_all.csv", respectively.

## Test code for imputation of missing distances 
```
Rscript test_imputation_H1.R 80
Rscript test_imputation_H2.R 80
```

Two partial trees are derived from "data_imputation/tree_1.txt".
These trees contain external nodes "data_imputation/idx1_L80.txt" and "data_imputation/idx2_L80.txt", respectively.
80 external nodes (+ outgroup node) are included in both of the partial trees.
The imputed distance matrix are saved in "out_imputation/optD_L80_H1.txt" and "out_imputation/optD_L80_H2.txt", respectively.



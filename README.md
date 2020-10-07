The script used for analyses in "Novel metric for hyperbolic phylogenetic tree embeddings"

## Reference

## Requirements

The following libraries are used.

* hydra
* MASS
* ape
* ROCR

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

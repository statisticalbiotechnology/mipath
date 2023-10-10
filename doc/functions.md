This package has two main functions: ``decompose_pathways`` and ``score_factors``. To get acquainted with their basic usage, please see the [MIPath vignette](https://github.com/statisticalbiotechnology/mipath/blob/main/MIPath_vignette.ipynb)

## decompose_pathways

```python
decompose_pathways(data, gene_sets_df, n_neighbors=25, min_genes=3, snn=True, snn_cutoff=0, knn_method = 'fast', n_jobs=-1)
```
* **data**: a pandas DataFrame containing the gene expression data. Rows are samples and columns are genes. Indexes are sample IDs and column names are gene IDs. This should be normalized and preprocessed before input.


* **gene_sets_df**: a pandas DataFrame containing the pathway annotation. The gene IDs must be the same as **data** column names. This should be obtained using either of the two *helper functions* below.

* **n_neighbors**: an integer parameter that determines the number of neighbors used to construct the graph. Default is 25.
The main function of this package is:

* **min_genes**: an integer parameter that is used to filter out small pathways. The size is measured in number of genes, and it considers only genes that are present both on the pathway annotation and the data. Default is 3.

* **snn**: a boolean, weather or not to wheith the nearest neighbors graph using the SNN similarity metric. Default is True and highly recomended.

* **snn_cutoff**: a float parameter, between 0 and 1, that is used to trim the SNN weighted graph, by removing edges with weight below the cutoff. Default is 0 (no cutoff), but it can be increased to speed up the algorithm, at the expense of better accuracy.

* **knn_method**: either 'fast' or 'exact', indicates which method will be used to build the NN graphs. Default is 'fast' and highly recomended.

* **n_jobs**: the number of parallel jobs that will be used in the 'fast' method for nearest neighbors. Default is -1 and indicates that all availiable cores will be used.


**Returns**: a pandas DataFrame, with same index as **data**, and with a column for each pathway in **gene_sets_df**.

## score_factors

```python
score_factors(decomposed_df, metadata, factors)
```
* **decomposed_df**: a pandas DataFrame, the result of the `decompose_pathways` function.

* **metadata**: a pandas DataFrame containing the phenotype annotations against witch the pathways will be scored. Indexes must be the same as **data** used in the `decompose_pathways` step.

* **factors**: a string or list of strings. Each string must be the same as a **metadata** column name. These are the variables that will be scored against the pathways.


**Returns**: a pandas DataFrame containing the results of the pathway analysis. Indexes are pathway IDs taken from **gene_sets_df** and column names are **factors**.


# Helper functions

## parse_gmt

```python
parse_gmt(gmt_path, gene_id_dict_path = None)
```

* **gmt_path**: the path to the '.gmt' file to be parsed.

* **gene_id_dict_path**: the path to a tab separated file. Default is **None**, but if a path to a file is provided, this file will be used to convert gene IDs. The file should be a *tab separated* file (tsv) containing two columns, the fist is a list of gene IDs in the original format, and the second the corresponding gene IDs in the desired format.

## get_reactome

```python
get_reactome(organism = 'HSA', gene_anot = 'Ensembl')
```

* **organism**: a string describing which organism to download the pathways for. Default is '**HSA**' for *Homo Sapiens*

* **gene_annot**: a string describing which gene annotation scheme to use. Accepts all Reactome schemes and defaults to '**Ensembl**'.

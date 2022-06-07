# Mipath documentation

## Package usage
### Main functions

The main function of this package is:

```python
mipath(data, metadata, gene_sets_df, factors, n_neighbors=25)
```

* **data**: a pandas DataFrame containing the gene expression data. Rows are samples and columns are genes. Indexes are sample IDs and column names are gene IDs. This should be normalized and preprocessed before input.
* **metadata**: a pandas DataFrame containing the phenotype annotations against witch the pathways will be scored. Indexes must be the same as **data**.
* **gene_sets_df**: a pandas DataFrame containing the pathway annotation. This should be obtained using the `parse_gmt` function in conjunction with `.gmt` file. The gene IDs must be the same as **data** column names. In case of reactome this can also be downloaded directly using the `get_reactome` helper function.
* **factors**: a string or list of strings. Each string must be the same as a **metadata** column name. These are the variables that will be scored against the pathways.
* **n_neighbors**: an integer parameter that determines the number of neighbors used to construct the graph. Default is 25.

This will return a pandas DataFrame containing the results of the pathway analysis. Indexes are pathway IDs taken from **gene_sets_df** and column names are **factors**.

The `mipath` function contains two steps, which can be run separately for convenience. Their separate functions are:

```python
decompose_pathways(data, gene_sets_df, n_neighbors=25)
```

which finds the sample modules for each pathway and returns a DataFrame containing the module assignment for each sample and pathway. Indexes are sample IDs and column names are pathway IDs. This is the most resource intensive part of the analysis.

The second step provides the final scores using:


```python
score_factors(decomposed_df, metadata, factors)
```

where **decomposed_df** is the result of the `decompose_pathways` function.

### Helper functions

```python
parse_gmt(gmt_path, gene_id_dict_path = None)
```

* **gmt_path**: the path to the `.gmt' file to be parsed.
* **gene_id_dict_path**: the path to a tab separated file. If not **None** this file will be used to convert gene IDs. It contains two columns, the fist is a list of gene IDs in the original format, and the second the corresponding gene IDs in the desired format.


```python
get_reactome(organism = 'HSA', gene_anot = 'Ensembl')
```

* **organism**: a string describing which organism to download the pathways for. Default is **HSA** for *Homo Sapiens*
* **gene_annot**: a string describing which gene annotation scheme to use. Accepts all Reactome schemes and defaults to **Ensembl**.

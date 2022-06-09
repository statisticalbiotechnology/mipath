import pandas as pd
import numpy as np

from sklearn.utils.multiclass import type_of_target
from sklearn.metrics.cluster import mutual_info_score

def significance_test(decomposed_df, metadata, factors, n_perms=1000):
    if not isinstance(factors, list):
        factors = [factors]
    # factors is a sting list containing the column headers
    result = pd.DataFrame(index = decomposed_df.columns, columns = factors)
    for label in factors:
        print("Significance test: ", label, ",", n_perms, "permutations")
        factor = metadata[label]
        if 'continuous' in type_of_target(factor.dropna()):
            print('Factor ' + label + ' is continuous')
            continue
        for pathway in decomposed_df.columns:
            temp_df = pd.DataFrame(columns=['pathway','factor'], index=decomposed_df.index)
            temp_df['pathway'] = decomposed_df[pathway]
            temp_df['factor'] = factor
            temp_df = temp_df.dropna(axis='index')
            pval = test(temp_df['pathway'], temp_df['factor'], n_perms)
            print(pathway + ": " + str(pval))
            result.loc[pathway, label] = pval
    return result

def test(X,Y, n_perms):
    mi = mutual_info_score(X,Y)
    null_mi = np.zeros(n_perms)
    for i in range(n_perms):
        X_perm = np.random.permutation(X)
        null_mi[i] = mutual_info_score(X_perm,Y)
    p = np.sum(null_mi >= mi) / n_perms
    return p
        


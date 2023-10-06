import pandas as pd
import numpy as np

from .pathways import get_reactome, get_pathway
from .decomposition import method_standard_decomposition
from .scoring import adjusted_contitional_mutual_information, adjusted_mutual_information

from joblib import Parallel, delayed, cpu_count

from sklearn.utils.multiclass import type_of_target

def decompose_pathways(data, gene_sets_df, n_neighbors=25, min_genes=3, snn=True, snn_cutoff=0, knn_method = 'fast', n_jobs=-1):
    print('Performing pathway decomposition')
    decomposed_df = pd.DataFrame(index = data.index, columns = gene_sets_df.index)
    i = 1
    for pathway in gene_sets_df.index:
        print('('+str(i)+'/'+str(len(gene_sets_df.index))+')','---',pathway,'---', len(gene_sets_df.loc[pathway, 'genes']), 'genes')
        i+=1
        try:
            data_path = get_pathway(data, gene_sets_df, pathway, min_genes=min_genes)
            data_path = data_path.dropna(axis=0)
        except ValueError:
            print('Not enough genes in pathway/data')
            print('')
            continue
        partition = method_standard_decomposition(data_path,
                n_neighbors=n_neighbors,
                snn_cutoff = snn_cutoff,
                knn_method=knn_method,
                metric = 'euclidean',
                leiden_iterations = -1,
                max_comm_size = 0,
                n_jobs=n_jobs)
        decomposed_df[pathway] = partition
        print('')
    return decomposed_df

def score_one_factor(decomposed_df, factor):
    #factor is a Series
    result = pd.Series(index = decomposed_df.columns)
    for pathway in decomposed_df.columns:
        if decomposed_df[pathway].isnull()[0] == False:
            temp_df = pd.DataFrame(columns=['pathway','factor'], index=decomposed_df.index)
            temp_df['pathway'] = decomposed_df[pathway]
            temp_df['factor'] = factor
            temp_df = temp_df.dropna(axis='index')
            score = adjusted_mutual_information(temp_df['pathway'], temp_df['factor'])
            result[pathway] = score
    return result

def score_one_factor_conditional(decomposed_df, factor, conditional):
    #factor and conditional are Series
    result = pd.Series(index = decomposed_df.columns)
    for pathway in decomposed_df.columns:
        if decomposed_df[pathway].isnull()[0] == False:
            temp_df = pd.DataFrame(columns=['pathway','factor', 'conditional'], index=decomposed_df.index)
            temp_df['pathway'] = decomposed_df[pathway]
            temp_df['factor'] = factor
            temp_df['conditional'] = conditional
            temp_df = temp_df.dropna(axis='index')
            score = adjusted_contitional_mutual_information(temp_df['pathway'], temp_df['factor'], temp_df['conditional'])
            result[pathway] = score
    return result


def score_proc(pathway_series, factor_series):
    if pathway_series.isnull().iloc[0] == True:
        return (pathway_series.name, np.nan)
    temp_df = pd.DataFrame(columns=['pathway','factor'], index=pathway_series.index)
    temp_df['pathway'] = pathway_series
    temp_df['factor'] = factor_series
    temp_df = temp_df.dropna(axis='index')
    score = adjusted_mutual_information(temp_df['pathway'], temp_df['factor'])
    return (pathway_series.name, score)

def score_one_factor_parallel(decomposed_df, factor):
    #factor is a Series
    ncores = cpu_count() - 1
    scores = Parallel(n_jobs=ncores)(delayed(score_proc)(pathway_series, factor) for pathway_id, pathway_series in decomposed_df.items())
    scores_df = pd.DataFrame(data = scores, columns=['pathway_id', factor.name]).set_index('pathway_id')[factor.name]
    return scores_df


def score_cond_proc(pathway_series, factor_series, conditional_series):
    if pathway_series.isnull()[0] == True:
        return (pathway_series.name, np.nan)
    temp_df = pd.DataFrame(columns=['pathway','factor', 'conditional'], index=pathway_series.index)
    temp_df['pathway'] = pathway_series
    temp_df['factor'] = factor_series
    temp_df['conditional'] = conditional_series
    temp_df = temp_df.dropna(axis='index')
    score = adjusted_contitional_mutual_information(temp_df['pathway'], temp_df['factor'], temp_df['conditional'])
    return (pathway_series.name, score)

def score_one_factor_conditional_parallel(decomposed_df, factor, conditional):
    #factor, conditional are pd.Series
    ncores = cpu_count() - 1
    scores = Parallel(n_jobs=ncores)(delayed(score_cond_proc)(pathway_series, factor, conditional) for pathway_id, pathway_series in decomposed_df.items())
    scores_df = pd.DataFrame(data = scores, columns=['pathway_id', factor.name]).set_index('pathway_id')[factor.name]
    return scores_df


def score_factors(decomposed_df, metadata, factors):
    # factors is a sting list containing the column headers
    if not isinstance(factors, list):
        factors = [factors]
    result = pd.DataFrame(index = decomposed_df.columns, columns = factors)
    for label in factors:
        print("Scoring", label)
        factor = metadata[label]
        if 'continuous' in type_of_target(factor.dropna()):
            print('Factor ' + label + ' is continuous')
            continue
        scores = score_one_factor_parallel(decomposed_df, factor)
        result[label] = scores
    return result

def score_factors_conditional(decomposed_df, metadata, factors, conditional):
    # factors is a sting list containing the column headers
    # conditional is pd.Series
    if not isinstance(factors, list):
        factors = [factors]
    result = pd.DataFrame(index = decomposed_df.columns, columns = factors)
    for label in factors:
        print("Scoring", label)
        factor = metadata[label]
        if 'continuous' in type_of_target(factor.dropna()):
            print('Factor ' + label + ' is continuous')
            continue
        scores = score_one_factor_conditional_parallel(decomposed_df, factor, conditional)
        result[label] = scores
    return result

def score_all_factors(decomposed_df, metadata):
    result = score_factors(decomposed_df, metadata, factors = metadata.columns)
    return result

def mipath(data, metadata, gene_sets_df, factors, n_neighbors=25):
    decomposed_df = decompose_pathways(data, gene_sets_df, n_neighbors)
    results = score_factors(decomposed_df, metadata, factors)
    return results


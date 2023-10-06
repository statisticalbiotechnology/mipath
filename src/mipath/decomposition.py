import pandas as pd
import numpy as np

import time
import scipy

from sklearn.neighbors import kneighbors_graph
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from itertools import repeat

import igraph as ig
import leidenalg as la
import pynndescent

def compute_nn(X, n_neighbors=20, metric='euclidean'):
    # Default parameters follow Seurat
    tic = time.time()
    knn_graph = kneighbors_graph(X, n_neighbors=n_neighbors, include_self=False, metric=metric, n_jobs=-1)
    print('Computed', n_neighbors, 'NN in', np.round(time.time()-tic,1), 's')
    return knn_graph

def compute_snn(knn_graph, cutoff=1/15):
    # Default parameters follow Seurat
    tic = time.time()
    n_neighbors = len(knn_graph[0].nonzero()[1])
    share_neighbors = knn_graph.dot(knn_graph.transpose())
    share_neighbors = share_neighbors/n_neighbors
    snn = share_neighbors >= cutoff
    print('SNN graph in',  np.round(time.time()-tic,1), 's')
    return snn

def snn_score(knn_graph, cutoff = 0):
    # Default parameters follow Seurat
    tic = time.time()
    n_neighbors = len(knn_graph[0].nonzero()[1])
    # Use adjecency power trick, uint8 for lower memory usage
    share_neighbors = knn_graph.astype(np.uint8).dot(knn_graph.astype(np.uint8).transpose())
    knn_scored = scipy.sparse.lil_matrix(knn_graph.shape)
    knn_scored[knn_graph.nonzero()] = share_neighbors[knn_graph.nonzero()] / n_neighbors
    knn_scored = knn_scored.tocsr()
    knn_scored.data[knn_scored.data < cutoff] = 0
    knn_scored.eliminate_zeros()
    print('SNN scores in',  np.round(time.time()-tic,1), 's')
    return knn_scored

def build_graph(matrix, directed=True):
    tic = time.time()
    # source, targets = matrix.nonzero()
    # g = ig.Graph(zip(source,targets), directed = directed)
    g = ig.Graph.Weighted_Adjacency(matrix)
    # print('iGraph in',  np.round(time.time()-tic,1), 's')
    return g

def knn_fast(X, n_neighbors=20, metric='euclidean', n_jobs=-1):
    # Default parameters follow Seurat
    tic = time.time()
    knn_search_index = pynndescent.NNDescent(X, n_neighbors=n_neighbors+1, metric=metric, n_jobs=n_jobs)
    knn_indices, knn_dists = knn_search_index.neighbor_graph
    knn_indices = knn_indices[:,1:] # Exclude self
    knn_graph = adjacency_matrix_representation(knn_indices, np.ones(knn_indices.shape)) 
    print('Computed', n_neighbors, 'NN in', np.round(time.time()-tic,1), 's')
    return knn_graph

def adjacency_matrix_representation(neighbor_indices, neighbor_distances):
    # Follows pyNNDescent function of the same name, but does not symmetrizes the output
    result = scipy.sparse.coo_matrix(
        (neighbor_indices.shape[0], neighbor_indices.shape[0]), dtype=np.float32
    )

    # Preserve any distance 0 points
    neighbor_distances[neighbor_distances == 0.0] = np.finfo(np.float32).eps

    result.row = np.repeat(
        np.arange(neighbor_indices.shape[0], dtype=np.int32), neighbor_indices.shape[1]
    )
    result.col = neighbor_indices.ravel()
    result.data = neighbor_distances.ravel()

    # Get rid of any -1 index entries
    result = result.tocsr()
    result.data[result.indices == -1] = 0.0
    result.eliminate_zeros()

    return result

def knn_exact(X, n_neighbors=20, metric='euclidean'):
    # Default parameters follow Seurat
    tic = time.time()
    knn_graph = kneighbors_graph(X, n_neighbors=n_neighbors, include_self=False, metric=metric, n_jobs=-2)
    print('Computed', n_neighbors, 'NN in', np.round(time.time()-tic,1), 's')
    return knn_graph

def find_partition(graph, n_iterations=-1, max_comm_size=0):
    tic = time.time()
    partition = la.find_partition(graph, la.ModularityVertexPartition, initial_membership=None, weights=graph.es['weight'], n_iterations=n_iterations, max_comm_size=max_comm_size)
    memberships = partition.membership
    n_part = len(np.unique(memberships))
    print('Found', n_part, 'partitions in', np.round(time.time()-tic,1), 's')
    return memberships

def method_standard_decomposition(in_df,
        n_neighbors=25,
        metric = 'euclidean',
        snn=True,
        snn_cutoff = 0,
        leiden_iterations = -1,
        max_comm_size = 0,
        knn_method = 'fast',
        n_jobs=-1):
    # Follows default parameters of Seurat
    if knn_method == 'fast':
        NN = knn_fast(in_df.values, n_neighbors=n_neighbors, metric=metric, n_jobs=n_jobs)
    elif knn_method == 'exact':
        NN = knn_exact(in_df.values, n_neighbors=n_neighbors, metric=metric)
    else:
        raise ValueError('Method ' + knn_method + ' not availiable')
    if snn:
        NN_scored = snn_score(NN, cutoff=snn_cutoff)
        G = build_graph(NN_scored, directed=True)
    else:
        G = build_graph(NN, directed=True)
    partition = find_partition(G, n_iterations=leiden_iterations, max_comm_size=max_comm_size)
    out_df = pd.Series(partition, index = in_df.index, dtype='category')
    return out_df


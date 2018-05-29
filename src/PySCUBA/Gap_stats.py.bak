#!/usr/bin/env python


# PySCUBA/src/PySCUBA/Gap_stats.py


# Author: Gregory Giecold for the GC Yuan Lab
# Affiliation: Harvard University
# Contact: g.giecold@gmail.com; ggiecold@jimmy.harvard.edu


from collections import defaultdict

import numpy as np
from sklearn.cluster import k_means
from sklearn.cluster import MiniBatchKMeans


__all__ = ['gap_stats']


def KMEANS(data, k):

    if data.shape[0] < 20000:
        centroids, cluster_IDs, _ = k_means(data, k, init = 'k-means++', precompute_distances = 'auto', n_init = 20, max_iter = 200)
    else:
        mbkm = MiniBatchKMeans(k, 'k-means++', max_iter = 100, batch_size = data.shape[0] / k, n_init = 20)
        mbkm.fit(data)
            
        centroids = mbkm.cluster_centers_
        cluster_IDs = mbkm.labels_

    return centroids, cluster_IDs


def box_corners(data):

    mins = np.amin(data, axis = 0)
    maxs = np.amax(data, axis = 0)

    return zip(mins, maxs)


def one_to_max(array_in):
    """Alter a vector of cluster labels to a dense mapping. 
       Given that this function is herein always called after passing 
       a vector to the function checkcl, one_to_max relies on the assumption 
       that cluster_run does not contain any NaN entries.

    Parameters
    ----------
    array_in : a list or one-dimensional array
        The list of cluster IDs to be processed.
    
    Returns
    -------
    result : one-dimensional array
        A massaged version of the input vector of cluster identities.
    """
    
    x = np.asanyarray(array_in)
    N_in = x.size
    array_in = np.reshape(x, N_in)    

    sorted_array = np.sort(array_in)
    sorting_indices = np.argsort(array_in)

    last = np.nan
    current_index = -1
    for i in xrange(N_in):
        if last != sorted_array[i] or np.isnan(last):
            last = sorted_array[i]
            current_index += 1

        sorted_array[i] = current_index

    result = np.empty(N_in, dtype = int)
    result[sorting_indices] = sorted_array

    return result


def W_k(data, centroids, cluster_IDs):

    cluster_IDs = one_to_max(cluster_IDs)

    N_clusters = int(np.amax(cluster_IDs) + 1)
    assert len(centroids) == N_clusters

    W_k = 0
    for i in xrange(N_clusters):
        samples_in_i = np.where(cluster_IDs == i)[0]
        if samples_in_i.size > 0:
            W_k += np.sum(np.linalg.norm(data[j] - centroids[i]) ** 2 for j in samples_in_i)

    return W_k


def gap_stats(data, min_k = 1, max_k = 10):

    B = 100

    assert isinstance(min_k, int) or type(min_k) is np.int_
    assert isinstance(max_k, int) or type(max_k) is np.int_
    assert (min_k > 0) and (max_k > 0)
    
    if min_k == max_k:
        return None, None, None
   
    k_list = np.arange(min_k, max_k + 1)

    min_maxs = box_corners(data)

    log_W_list = []
    E_log_W_list = []
    s_k_list = []
    
    for k in k_list:
        centroids, cluster_IDs = KMEANS(data, k)

        log_W_list.append(np.log(W_k(data, centroids, cluster_IDs))) 

        log_W_k_b_list = np.zeros(B, dtype = float)
        for b in xrange(B):
            uniform_data = np.zeros((data.shape[0], data.shape[1]), dtype = float)
            for i in xrange(data.shape[1]):
                v = np.random.uniform(low = min_maxs[i][0], high = min_maxs[i][1], size = data.shape[0])
                uniform_data[:, i] = v

            centroids, cluster_IDs = KMEANS(uniform_data, k)

            log_W_k_b_list[b] = np.log(W_k(uniform_data, centroids, cluster_IDs))            

        E_log_W_list.append((np.sum(log_W_k_b_list) + 0.0) / B)

        s_k = np.sum((log_W_k_b_list - E_log_W_list[-1]) ** 2) / B
        s_k = np.sqrt(s_k)
        s_k_list.append(s_k)

    log_W_list = np.asarray(log_W_list, dtype = float)     
    E_log_W_list = np.asarray(E_log_W_list, dtype = float)

    s_k_list = np.asarray(s_k_list, dtype = float)
    s_k_list *= np.sqrt(1+ 1/(B + 0.0))

    return log_W_list, E_log_W_list, s_k_list


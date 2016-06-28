#!/usr/bin/env python


# PySCUBA/src/PySCUBA/SCUBA_core.py


# Author: Gregory Giecold for the GC Yuan Lab
# Affiliation: Harvard University
# Contact: g.giecold@gmail.com; ggiecold@jimmy.harvard.edu


from collections import defaultdict
from math import exp, floor, sqrt
from os import getcwd, path
import sys
import warnings

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import scipy.optimize
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances_argmin
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.preprocessing import StandardScaler

from .Gap_stats import gap_stats


__all__ = ['bifurcation_analysis', 'bifurcation_direction', 'initialize_tree', 
           'PCA_analysis', 'refine_tree']


SCALING_FACTOR = float(10)


class Eval:
    """This class conveniently enables the evaluation inside strings of expressions
       involving local and global variables.
    """
    
    def __init__(self, globals = None, locals = None):
        self.globals_ = globals or {}
        self.locals_ = locals

    def __getitem__(self, key):
        if self.locals_ is None:
            self.locals_ = sys._getframe(1).f_locals
            
        if len(self.globals_) == 0:
            self.globals_ = sys._getframe(1).f_globals

        key = key % self
        
        return eval(key, self.globals_, self.locals_)


def PCA_analysis(data, mode, cell_stages = None):
    """Principal Component Analysis.
    """

    assert mode in {'pca', 'pca2'}
    
    mean_shifter = StandardScaler(with_std = False)
    
    if mode == 'pca':
        pca = PCA(min(data.shape))
        projected_data = pca.fit_transform(data)
        projected_data = pca.fit_transform(mean_shifter.fit_transform(data))
        components = pca.components_
    else:
        assert isinstance(cell_stages, np.ndarray)
        
        idx = np.where(cell_stages == np.max(cell_stages))[0]
        
        pca = PCA(min(idx.size, data.shape[1]))
        pca.fit(mean_shifter.fit_transform(data[idx]))
        components = pca.components_
        projected_data = np.dot(data, components.T)
    
    return components, projected_data
    
    
def custom_formatwarning(msg, *a):
    """Display a custom message and ignore every other warning.
    """
    
    return str(msg) + '\n'
    
    
def initialize_tree(data, cell_stages, rigorous_gap_stats = False):
    """Come up with an initial assessment of the cellular hierarchy, using a series of
       k-means clustering.
    """

    warnings.formatwarning = custom_formatwarning

    min_split = 15                 # Lower threshold on the number of cells in a cluster 
                                   # for this cluster to be split.
    min_percentage_split = 0.25     # Minimum fraction of cells in the smaller cluster
                                   # during a bifurcation.
    
    N_cells, N_features = data.shape
    
    stages = np.unique(cell_stages)
    N_stages = stages.size
    
    cluster_indices = np.full(N_cells, -1, dtype = int)
    parent_clusters = defaultdict(list)
    cluster_IDs_per_stage = defaultdict(list)
    
    # We're now about to compute a gap statistics estimation of the optimal number
    # of clusters at the initial (pseudo-)time stage:
    
    condition = cell_stages == stages[0]
    
    X = np.compress(condition, data, axis = 0)
    
    if rigorous_gap_stats and (X.shape[0] < min_split):
        N_best_clusters = 1
    elif rigorous_gap_stats:
        N_best_clusters = 1
        N_best_clusters = 1
        k_max = min(2, X.shape[0])
    
        log_W, E_log_W, s_k = gap_stats(X, 1, k_max)
        if log_W is not None:
            gaps = E_log_W - log_W
            for l in xrange(gaps.size - 1):
                if gaps[l] >= gaps[l + 1] - s_k[l + 1]:
                    N_best_clusters = l + 1
                    break
    else:
        N_samples = 1000
        k_max = min(2, X.shape[0])
        
        R = np.random.uniform(size = (N_samples, N_features))
        W = np.diagflat(np.std(data, axis = 0))
        R = np.dot(R, W) / sqrt(1 / float(12))
   
        W_k0 = np.empty(k_max, dtype = float)
        for i in xrange(k_max):
            kmeans = KMeans(i + 1, n_init = 50)
            W_k0[i] = kmeans.fit(R).inertia_ / float(N_samples)
    
        W_k = np.empty(k_max, dtype = float)
        for i in xrange(k_max):
            kmeans = KMeans(i + 1, n_init = 50)
            W_k[i] = kmeans.fit(X).inertia_ / float(X.shape[0])
        
        N_best_clusters = np.argmax(W_k0[:k_max] - W_k[:]) + 1
        if X.shape[0] < min_split:
            N_best_clusters = 1
         
    # Done with the gap statistics for the first (pseudo-)time stage
    
    cluster_tally = N_best_clusters
    
    kmeans = KMeans(N_best_clusters, n_init = 50)
    kmeans.fit(X)
    
    cluster_indices[condition] = kmeans.labels_
    centroid_coordinates = kmeans.cluster_centers_
    parent_clusters[0].extend([-1] * N_best_clusters)
    cluster_IDs_per_stage[0].extend([k for k in xrange(N_best_clusters)])
    
    for stage_idx in xrange(1, N_stages):
        condition = cell_stages == stages[stage_idx]
    
        X = np.compress(condition, data, axis = 0)
        
        # Map the cells at this stage to the nearest clusters from the previous stage:
        previous_clusters = cluster_IDs_per_stage[stage_idx - 1]
        N_previous_clusters = len(previous_clusters)
     
        previous_centroid_coordinates = centroid_coordinates[-N_previous_clusters:]
        nearest_previous_cluster_indices = pairwise_distances_argmin(X,
                                          previous_centroid_coordinates)
        
        # The following will now test for a bifurcation at this (pseudo-)time stage:
        for j in xrange(N_previous_clusters):
            idx = np.where(nearest_previous_cluster_indices == j)[0]
            
            if idx.size == 0:
                msg = ' '.join(["Empty cluster mapping encountered", 
                                "at (pseudo-)time {0}.\n".format(stage_idx + 1)])
                warnings.warn(msg)
                continue
                
            XX = X[idx]
            
            if rigorous_gap_stats:
                N_best_clusters = 1
                k_max = min(2, XX.shape[0])
            
                log_W, E_log_W, s_k = gap_stats(XX, 1, k_max)
                if log_W is not None:
                    gaps = E_log_W - log_W
                    if gaps.size > 1 and gaps[0] >= gaps[1] - s_k[1]:
                        N_best_clusters = 2
                        kmeans = KMeans(2, n_init = 50)
                        cluster_idx = kmeans.fit(XX).labels_
            else:
                W_k = []
                for i in {1, 2}:
                    if XX.shape[0] >= i:
                        kmeans = KMeans(i, n_init = 50)
                        W_k.append(kmeans.fit(XX).inertia_ / float(XX.shape[0]))
                        
                N_best_clusters = np.argmax(W_k0[:len(W_k)] - W_k[:]) + 1
                kmeans = KMeans(N_best_clusters, n_init = 50)
                cluster_idx = kmeans.fit(XX).labels_
            
            # Do not split if there are not enough cells 
            if (idx.size < min_split) or (N_best_clusters == 1) or (min(np.bincount(cluster_idx)) / float(XX.shape[0]) < min_percentage_split):
                cluster_indices[np.where(condition == True)[0][idx]] = cluster_tally
                centroid_coordinates = np.vstack((centroid_coordinates, np.mean(XX, 0)))
                parent_clusters[stage_idx].append(cluster_IDs_per_stage[stage_idx-1][j])
                cluster_IDs_per_stage[stage_idx].append(cluster_tally)
                
                cluster_tally += 1
            else:
                cluster_indices[np.where(condition == True)[0][idx]] = cluster_tally + cluster_idx
                centroid_coordinates = np.vstack((centroid_coordinates,
                                               kmeans.cluster_centers_))
                parent_clusters[stage_idx].extend([cluster_IDs_per_stage[stage_idx-1][j]] * 2)
                cluster_IDs_per_stage[stage_idx].extend([cluster_tally, cluster_tally + 1])
                
                cluster_tally += 2
                
    assert np.all(cluster_indices != -1)
                
    return centroid_coordinates, cluster_indices, parent_clusters


def refine_tree(data, centroid_coordinates, cluster_indices, parent_clusters, cell_stages, 
                output_directory = None):
    """Refine the tree structure initially determined by the procedure 'initial_tree' above
       by maximizing the penalized likelihood function [Eq. (1) in Marco E., Karp R.L., Guo G.,
       Robson P., Hart A.H., Trippa L. and Yuan, G.-C., 'Bifurcation Analysis of Single-Cell Gene
       Expression Data Reveals Epigenetic Landscape.' In: Proc. Natl. Acad. Sci. USA, 2014 Dec. 30,
       111 (52)]
    """

    warnings.formatwarning = custom_formatwarning

    if output_directory is None:
        output_directory = getcwd()
        
    Lambda = 0.5    # Diffusion constant for cluster drift.

    N_cells, N_features = data.shape

    stages = np.unique(cell_stages)
    N_stages = stages.size
    
    N_clusters = centroid_coordinates.shape[0]
    
    N_clusters_per_stage = np.array([len(v) for v in parent_clusters.values()], dtype = int)
    cumulative_N_clusters = np.append([0], np.cumsum(N_clusters_per_stage))
    
    # Assume that the variances of each features are identical and the same across 
    # (pseudo-)time stages and clusters:
    Sigma = np.ones((N_clusters, N_features), dtype = float)
    
    # Starting the procedure that will update the tree structure based
    # the computation of a penalized likelihood:
    
    convergence_threshold = 10 ** (-8)    # Criterion for convergence.
        
    max_iter = 100    # Maximum number of iterations
    iter = 0
    while iter < max_iter:
        outdated_cluster_indices = np.copy(cluster_indices)
        outdated_centroid_coordinates = np.copy(centroid_coordinates)
        
        outdated_parent_clusters = defaultdict(list)
        for k, v in parent_clusters.items():
            outdated_parent_clusters[k] = v
        
        for stage_idx in xrange(N_stages):
            beg = cumulative_N_clusters[stage_idx]
            end = cumulative_N_clusters[stage_idx + 1]
            
            centroid_coordinates_at_stage = centroid_coordinates[beg:end]
            Sigma_at_stage = Sigma[beg:end]
                
            cells_at_stage = np.where(cell_stages == stages[stage_idx])[0]
            
            if end - beg > 1:
                dispersion = np.empty((end - beg, cells_at_stage.size), dtype = float)
                for i in xrange(end - beg):
                    tmp = data[cells_at_stage] - centroid_coordinates[beg:end][i]
                    tmp = np.divide(tmp, Sigma_at_stage[i].astype(float))
                    tmp = tmp ** 2
                    dispersion[i] = np.sqrt(np.mean(tmp, axis = 1))
                
                cluster_idx = np.argmin(dispersion, axis = 0)
                cluster_indices[cells_at_stage] = beg + cluster_idx
                
                empty_clusters = np.where(np.bincount(cluster_idx) == 0)[0]
                if empty_clusters.size > 0:
                    msg = ' '.join(["WARNING: PySCUBA: refine_tree: empty cluster",
                        "identified at (pseudo-)time level {0}.\n".format(stage_idx + 1)])
                    warnings.warn(msg)
                    
                    used_samples = set()
                    for cluster in empty_clusters:
                        tmp = dispersion[cluster].reshape(dispersion.shape[1])
                        if len(used_samples):
                            tmp[list[used_samples]] = np.nan
                        
                        sample_idx = np.nanargmin(tmp)
                        used_samples.add(sample_idx)
                            
                        cluster_indices[cells_at_stage[sample_idx]] = beg + cluster
            # Updating the centroids:
            for i, cluster in enumerate(cumulative_N_clusters[stage_idx:stage_idx + 1]):
                samples_in_cluster = np.where(cluster_indices == cluster)[0]
                
                tmp = np.sum(data[samples_in_cluster], axis = 0)
                
                parent_idx = parent_clusters[stage_idx][i]
                if parent_idx == -1:
                    centroid_coordinates[cluster] = tmp / float(tmp.size)
                else:
                    tmp = 2 * Lambda * outdated_centroid_coordinates[parent_idx] + tmp
                    centroid_coordinates[cluster] = tmp / float(2 * Lambda + tmp.size)
                    
            # Updating the parental clusters:
            if stage_idx > 0:
                beg_previous = cumulative_N_clusters[stage_idx - 1]
                end_previous = cumulative_N_clusters[stage_idx]
                previous_centroid_coordinates = centroid_coordinates[beg_previous:end_previous]
                
                dispersion = euclidean_distances(previous_centroid_coordinates, centroid_coordinates_at_stage, squared = True)
                
                reshuffled_parent_clusters = np.argmin(dispersion, axis = 0)
                tmp = np.bincount(reshuffled_parent_clusters)
                if min(tmp) > 0 and max(tmp) < 3:
                    parent_clusters[stage_idx] = list(reshuffled_parent_clusters + beg_previous)
                else:
                    msg = ' '.join(["WARNING: PySCUBA: refine_tree: Binary tree",
                        "assumption violated at stage {0}.".format(stage_idx), 
                        "Using the unaltered configuration instead.\n"])
                    warnings.warn(msg)
                
        # Check if convergence has been reached yet:
        delta_centroids = np.sum(np.fabs(centroid_coordinates - outdated_centroid_coordinates))
        delta_indices = np.sum(np.fabs(cluster_indices - outdated_cluster_indices))
        
        delta_parents = 0
        for v1, v2 in zip(outdated_parent_clusters.values(), parent_clusters.values()):
            delta_parents += np.sum(np.fabs(np.subtract(v1, v2)))
        
        if iter > 20 and delta_centroids < convergence_threshold and delta_indices == 0 and delta_parents == 0:
            print("INFO: PySCUBA: refine_tree: iterations successfully converged "
                  "to an end.\n")
            break

        iter += 1
       
    
    # Writing to file the structure of the final SCUBA tree:
    final_tree_path = path.join(output_directory, 'final_tree.tsv')
        
    stages = np.empty(0, dtype = int)
    parents = np.empty(0, dtype = int)
    for k, v in parent_clusters.items():
        stages = np.concatenate((stages, [k] * len(v)))
        parents = np.append(parents, np.add(1, v))
        
    stages = stages.reshape(stages.size, 1)
    parents = parents.reshape(parents.size, 1)
        
    with open(final_tree_path, 'w') as f:
        f.write('\t'.join(['Cluster ID', 'Stage', 'Parent cluster\t']))
        f.write('\t'.join(['MU {0}'.format(k + 1) for k in xrange(N_features)]))
        f.write('\n')
        
    final_tree = np.genfromtxt(final_tree_path, delimiter = '\t', dtype = str)
    bundle = np.hstack((np.arange(1, N_clusters + 1).reshape(N_clusters, 1), stages, 
                       parents, centroid_coordinates))
    final_tree = np.vstack((final_tree, bundle))
    
    with open(final_tree_path, 'w') as f:
        np.savetxt(f, final_tree, fmt = '%s', delimiter = '\t')
        
    with open(path.join(output_directory, 'cluster_IDs.csv'), 'w') as f:
        np.savetxt(f, cluster_indices, fmt = '%d', delimiter = ',')

    return centroid_coordinates, cluster_indices, parent_clusters


def bifurcation_direction(data, cell_IDs, markers, parent_clusters, centroid_coordinates,
                          output_directory = None, weights = None):
    """Infer the direction associated with each bifurcation and project the data onto 
       the bifurcation directions.
    """

    warnings.formatwarning = custom_formatwarning

    if output_directory is None:
        output_directory = getcwd()
        
    N_clusters_per_stage = np.array([len(v) for v in parent_clusters.values()], dtype = int)
    cumulative_N_clusters = np.append([0], np.cumsum(N_clusters_per_stage))
    
    bifurcation_info = defaultdict(list)
    bifurcation_directions = np.empty((markers.size, 0), dtype = float)
    bifurcation_projections = np.empty((cell_IDs.size, 0), dtype = float)
    
    for stage_idx, parents in parent_clusters.items():
        if stage_idx > 0:
            tmp1 = np.bincount(parents)
            splitting_parents = np.where(tmp1 > 1)[0]
            for parent in splitting_parents:
                grandpa = parent_clusters[stage_idx - 1][parent - cumulative_N_clusters[stage_idx - 1]]
                children = list(np.where(parents == parent)[0] + cumulative_N_clusters[stage_idx])
                
                bifurcation_info[stage_idx].append([grandpa, parent] + children)
                
                tmp2 = centroid_coordinates[children[0]] - centroid_coordinates[children[1]]
                tmp2 /= float(np.dot(tmp2, tmp2))
                tmp2 *= SCALING_FACTOR
                tmp2 = tmp2.reshape(tmp2.size, 1)
                
                bifurcation_projections = np.hstack((bifurcation_projections, np.dot(data, tmp2)))
                
                if weights is None:
                    bifurcation_directions = np.hstack((bifurcation_directions,
                                                        tmp2))
                else:
                    bifurcation_directions = np.hstack((bifurcation_directions, 
                                                        np.dot(weights.T, tmp2)))
                    
    header = ['Stage']
    for k, v in bifurcation_info.items():
        header.extend([str(k)] * len(v))
    
    if len(header) == 1:
        warnings.warn('No bifurcation found!\n')
        return (None, None, None)
    
    bifurcation_data_path = path.join(output_directory, 'bifurcation_direction.tsv')
    with open(bifurcation_data_path, 'w') as f:
        f.write('\t'.join(header))
        f.write('\n')
        
    bifurcation_data = np.genfromtxt(bifurcation_data_path, delimiter = '\t', dtype = str).reshape(1, len(header))
    tmp = np.hstack((markers.reshape(markers.size, 1), bifurcation_directions))
    bifurcation_data = np.vstack((bifurcation_data, tmp))
    
    with open(bifurcation_data_path, 'w') as f:
        np.savetxt(f, bifurcation_data, fmt = '%s', delimiter = '\t')

    projection_data_path = path.join(output_directory, 'projection_all_data.tsv')
    with open(projection_data_path, 'w') as f:
        f.write('\t'.join(header))
        f.write('\n')
        
    projection_data = np.genfromtxt(projection_data_path, delimiter = '\t', dtype = str).reshape(1, len(header))
    tmp = np.hstack((cell_IDs.reshape(cell_IDs.size, 1), bifurcation_projections))
    projection_data = np.vstack((projection_data, tmp))
    
    with open(projection_data_path, 'w') as f:
        np.savetxt(f, projection_data, fmt = '%s', delimiter = '\t')
                                
    return bifurcation_info, bifurcation_directions, bifurcation_projections
                                   
      
def bifurcation_analysis(cluster_indices, bifurcation_info, bifurcation_directions,
                         bifurcation_projections, output_directory = None,
                         normalize_likelihood_level_cell_counts = True):
    """Infer the dynamical changes of gene expression patterns along each bifurcation direction
       by fitting a Fokker-Planck equation. Proceed by considering the gene expression data
       for several cell stages projected along a bifurcation axis. Perform a maximum likelihood
       estimation of the parameters of a Fokker-Planck equation with quartic potential 
       from the data at each such stage.
       
       Parameters
       ----------
       normalize_likelihood_level_cell_counts : bool, optional (default = True)
           Option to compute the log-likelihood with the average number of cells per cell stage
    """   
    
    warnings.formatwarning = custom_formatwarning
    
    assert bifurcation_directions is not None
    
    if output_directory is None:
        output_directory = getcwd()
        
    parameters_path = path.join(output_directory, 'Fokker_Planck_parameters.tsv')
    with open(parameters_path, 'w') as f:
        f.write('Stage:\tSigma:\tb:\ta_before:\ta_middle:\ta_after:\n')
    
    min_cells_per_bifurcation = 15
    
    data_per_split = dict()
    parameters_per_split = dict()
    counter = 0
    for stage, vectors in sorted(bifurcation_info.items(), key = lambda tpl: tpl[0]):
        for v in vectors:
            projected_data = bifurcation_projections[:, counter]
        
            data_grandpa = np.compress(cluster_indices == v[0], projected_data)
            data_grandpa -= np.mean(data_grandpa) 
            
            data_parent = np.compress(cluster_indices == v[1], projected_data)
            data_parent -= np.mean(data_parent)
            
            data_left_child = np.compress(cluster_indices == v[2], projected_data)
            data_right_child = np.compress(cluster_indices == v[3], projected_data)
            data_center = 0.5 * (np.mean(data_left_child) + np.mean(data_right_child))
            data_left_child -= data_center
            data_right_child -= data_center
            
            data = (data_grandpa, data_parent, np.append(data_left_child, data_right_child))
            
            ordinals = get_ordinals(bifurcation_directions.shape[1])
            
            if (2 * data[0].size > min_cells_per_bifurcation) and (2 * data[1].size > min_cells_per_bifurcation) and (data[2].size > min_cells_per_bifurcation):
                seed_a_parameter = (0.5 * (np.mean(data_right_child) - np.mean(data_left_child)))**2
                seed_parameters = (np.var(data[-1]), 0, -2 * seed_a_parameter, 
                                   -seed_a_parameter, seed_a_parameter)
                
                sys.stdout.write("\nEstimating the parameters of a Fokker-Planck model "
                                 "by maximizing the underlying likelihood functions associated "
                                 "with the {0} bifurcation.".format(ordinals[counter]))
                sys.stdout.write('\n...')
                sys.stdout.flush()
                                 
                parameters = MLE(data, seed_parameters, normalize_likelihood_level_cell_counts)
                
                sys.stdout.write('\nMLE completed for this bifurcation.\n')
                sys.stdout.flush()
                
                # Writing the fitted Fokker-Planck parameters to disk:
                parameters_string = '\t'.join([str(round(p, 4)) for p in parameters])
                with open(parameters_path, 'a') as f:
                    f.write("{stage}\t".format(**locals()) + parameters_string + '\n')
                    
                plot_potentials(counter, data, parameters, output_directory)
                
                data_per_split[counter] = data
                parameters_per_split[counter] = parameters
            else: 
                msg = ' '.join(["Not enough cells to fit the bifurcation potential", 
                    "of the {0} bifurcation at stage {1}.".format(ordinals[counter],
                    stage), "At least {0} cells are required.\n".format(
                    min_cells_per_bifurcation)])
                warnings.warn(msg)
            
            counter += 1
            
    return data_per_split, parameters_per_split


def normalize(sigma, a, b, limit):

    f = lambda y: exp((y**2 * (a - 0.5 * y**2) + 2 * b * y) / sigma**2)
    
    limit = abs(limit)
    normalization_factor = 1.0 / scipy.integrate.quad(f, -limit, limit)[0]
    
    return round(normalization_factor, 7)

    
def steady_state_distribution(x, sigma, a, b, normalization_factor):
    
    f = lambda y: exp((y**2 * (a - 0.5 * y**2) + 2 * b * y) / sigma**2)
    
    return round(f(x) * normalization_factor + (10 ** 10) * sys.float_info.epsilon, 7)
    

def multistate_distribution(data, parameters, limit, 
                            normalize_likelihood_level_cell_counts = True):

    data_grandpa, data_parent, data_children = data
    sigma, b, a_grandpa, a_parent, a_children = parameters
    
    normalization_factor = normalize(sigma, a_grandpa, b, limit)
    grandpa_dist = [steady_state_distribution(x, sigma, a_grandpa, b, normalization_factor) for x in data_grandpa]
    
    normalization_factor = normalize(sigma, a_parent, b, limit)
    parent_dist = [steady_state_distribution(x, sigma, a_parent, b, normalization_factor) for x in data_parent]
    
    normalization_factor = normalize(sigma, a_children, b, limit)
    children_dist = [steady_state_distribution(x, sigma, a_children, b, normalization_factor) for x in data_children]
    
    grandpa_dist = np.array(grandpa_dist, dtype = float)
    parent_dist = np.array(parent_dist, dtype = float)
    children_dist = np.array(children_dist, dtype = float)
    
    if normalize_likelihood_level_cell_counts:
        grandpa_dist = np.divide(grandpa_dist, float(data_grandpa.size))
        parent_dist = np.divide(parent_dist, float(data_parent.size))
        children_dist = np.divide(children_dist, float(data_children.size))
        
    return grandpa_dist, parent_dist, children_dist


def MLE(data, seed_parameters, normalize_likelihood_level_cell_counts = True):
    
    limit = int(floor(1.5 * np.max(np.abs(data[-1]))))
    
    def negative_log_likelihood(parameters):
        distributions = multistate_distribution(data, parameters, limit,
                                  normalize_likelihood_level_cell_counts)
                                       
        NLL = 0
        for dist in distributions:
            NLL -= np.sum(np.log(dist))
            
        return NLL
        
    bounds = ((0.01, None), (None, None), (None, 0), (None, 0), (0, None))
        
    results = scipy.optimize.minimize(negative_log_likelihood, seed_parameters, 
                                      method = 'SLSQP', bounds = bounds, 
                                      options = {'maxiter': 10000})
                              
    return results.x


def potential(x, parameters):

    sigma, b, a = parameters

    return -0.5 * (x**2 * (a - 0.5 * x**2) + 2 * b * x) / sigma**2
    
    
def get_ordinals(N):

    ordinals = ['1st', '2nd', '3rd']
    ordinals.extend(['{0}th'.format(i) for i in xrange(4, N + 1)])
    
    return ordinals
    
    
def plot_potentials(count, projected_data, parameters, output_directory = None):

    # Smoothing parameter for the histograms generated from the probability density functions
    # inferred one node before a bifurcation, at the bifurcation and right after it.
    N_bins = 40
    
    if output_directory is None:
        output_directory = getcwd()
    
    ordinals = get_ordinals(count + 1)
    
    limit = int(floor(1.5 * np.max(np.abs(projected_data[-1]))))
    grid = np.linspace(-limit, limit, num = 5000)
    
    fig = plt.figure(1, (11, 9))
    gs = gridspec.GridSpec(3, 2)

    sigma, b, a_grandpa, a_parent, a_children = parameters
    sigma = round(sigma, 3)
    b = round(b, 3)
    a_parameters = [a_grandpa, a_parent, a_children]
        
    for i, data in enumerate(projected_data):
        ax = plt.subplot(gs[i, 1])
        
        hist, bin_edges = np.histogram(data, bins= np.linspace(-limit, limit, N_bins))
        width = 0.7 * (bin_edges[1] - bin_edges[0])
        center = (bin_edges[:-1] + bin_edges[1:]) / 2
        ax.bar(center, hist, align = 'center', width = width, color = 'b')
            
        a = a_parameters[i]
        normalization_factor =  normalize(sigma, a, b, limit)
        fitted_distribution = np.array([steady_state_distribution(x, sigma, a, b, normalization_factor) for x in grid], dtype = float)
        fitted_distribution *= hist.max() / float(fitted_distribution.max())
        fitted_distribution += 1.2 
        ax.plot(grid, fitted_distribution, color = 'r')
            
        x_lim = int(floor(0.8 * limit))    
        ax.set_xticks(np.linspace(-x_lim, x_lim, 7))
        ax.set_xticklabels([str(round(k,1)) for k in np.linspace(-x_lim / SCALING_FACTOR, x_lim / SCALING_FACTOR, 7)])
        ax.set_yticks(np.linspace(0, hist.max() + 1, 7))
            
        ax.set_xlim([-limit, limit])
        ax.set_ylim([0, hist.max() + 2])
        
        if 2 == i:    
            ax.set_xlabel('Scores on bifurcation axes')
        ax.set_ylabel('Counts')
        if 0 == i:
            ax.set_title('Estimate of PDF and histogram')

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(9)

        fig.add_subplot(ax)
       
    for i, a in enumerate(a_parameters):
        a = round(a, 3)
            
        p = np.array([potential(x, (sigma, b, a)) for x in grid], dtype = float)
            
        ax = plt.subplot(gs[i, 0])
        ax.plot(grid, p, color = 'm')
        
        x_lim = int(floor(0.8 * limit))    
        ax.set_xticks(np.linspace(-x_lim, x_lim, 7))
        ax.set_xticklabels([str(round(k,1)) for k in np.linspace(-x_lim / SCALING_FACTOR, x_lim / SCALING_FACTOR, 7)])
        ax.set_yticks(np.linspace(p.min() - 1.9, int(floor(np.median(p))) + 2.9, 7))
            
        ax.set_xlim([-limit, limit])
        ax.set_ylim([p.min() - 2, int(floor(np.median(p))) + 3])
        
        if 2 == i:    
            ax.set_xlabel('Scores on bifurcation axes')
        ax.set_ylabel('Potential')
        ax.set_title(r"$\sigma$ = %(sigma / SCALING_FACTOR).3f, b = %(b / SCALING_FACTOR).3f, a = %(a / SCALING_FACTOR).3f" % Eval())
            
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(9)
                
        fig.add_subplot(ax)
        
    plt.savefig(path.join(output_directory, 'fitted_potentials_for_{}_bifurcation.pdf'.format(ordinals[count])))
    plt.close(fig)


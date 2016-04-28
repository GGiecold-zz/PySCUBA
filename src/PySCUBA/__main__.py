#/usr/env/bin python


# PySCUBA/src/PySCUBA/__main__.py


# Author: Gregory Giecold for the GC Yuan Lab
# Affiliation: Harvard University
# Contact: g.giecold@gmail.com; ggiecold@jimmy.harvard.edu


import optparse
from os import getcwd, path
import sys

import numpy as np
from sklearn.preprocessing import StandardScaler

from Gap_stats import gap_stats
from Preprocessing import cytometry_preprocess, PCR_preprocess, RNASeq_preprocess
import .SCUBA_core as SCUBA


def plot_tree(cluster_indices, parent_clusters, output_directory = None):
    """Display a bifurcation tree.
    """

    if output_directory is None:
        output_directory = getcwd()
        
    vertex_sizes = np.bincount(cluster_indices)
    N_vertices = vertex_sizes.size
    vertex_sizes = np.divide(vertex_sizes, float(np.sum(vertex_sizes)))
    vertex_sizes *= 100 * N_vertices
    vertex_sizes += 40 + (N_vertices / 3)

    tree = igraph.Graph()
    tree.add_vertices(N_vertices)
    
    cluster_tally = 0
    for k, v in parent_clusters.items():
        if k > 0:
            tree.add_edges(zip(v, xrange(cluster_tally, cluster_tally + len(v))))
        cluster_tally += len(v)
        
    tree.vs['label'] = xrange(N_vertices)
    layout = tree.layout('fr')
    name = path.join(output_directory, 'SCUBA_tree.pdf')
    igraph.plot(tree, name, bbox = (200 * N_vertices, 200 * N_vertices), margin = 250,
                layout = layout, edge_width = [7] * (N_vertices - 1), vertex_label_dist = 0, 
                vertex_label_size = 30, vertex_size = vertex_sizes.tolist())
    

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
    array_in = x.reshape(N_in)    

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
       
              
def main():
        
    opts, file_path = parser.parse_args()

    preprocessing_fcts = [cytometry_preprocess, PCR_preprocess, RNASeq_preprocess]
    data_type = {'cytometry': 0, 'PCR': 1, 'RNASeq': 2}
    cell_IDs, data, markers, cell_stages, data_tag, output_directory = preprocessing_fcts[data_type[opts.data_type]](file_path, opts.log_mode, opts.pseudotime_mode, opts.pcv_method, opts.anchor_gene, opts.exclude_marker_names)
    
    cell_stages = 1 + one_to_max(cell_stages)
    data = StandardScaler(with_std = False).fit_transform(data)
    
    cluster_mode = opts.cluster_mode
    
    if cluster_mode in {'pca', 'pca2'}:    
        PCA_components, data = SCUBA.PCA_analysis(data, cluster_mode, 
                    cell_stages if (cluster_mode == 'pca2') else None)
    
    centroid_coords, cluster_indices, parent_clusters = SCUBA.initialize_tree(data, cell_stages)
    centroid_coords, cluster_indices, parent_clusters = SCUBA.refine_tree(data, 
                             centroid_coords, cluster_indices, parent_clusters, 
                             cell_stages, output_directory)
    
    try:
        import igraph
        plot_tree(cluster_indices, parent_clusters, output_directory)
    except ImportError, e:
        pass
        
    if cluster_mode in {'pca', 'pca2'}:
        weights = PCA_components
    else:
        weights = None
    
    bifurcation_info, bifurcation_axes, bifurcation_projections = SCUBA.bifurcation_direction(data,
                    cell_IDs, markers, parent_clusters, centroid_coords, output_directory, weights)
                
    if bifurcation_info:
        data_per_split, parameters_per_split = SCUBA.bifurcation_analysis(cluster_indices, 
                              bifurcation_info, bifurcation_axes, bifurcation_projections,
                              output_directory)
        SCUBA.reduction_simulations(data_per_split, parameters_per_split)
    

if __name__ == '__main__':

    main()
    
    

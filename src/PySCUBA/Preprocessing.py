#!/usr/bin/env python


# PySCUBA/src/PySCUBA/Preprocessing.py


# Author: Gregory Giecold for the GC Yuan Lab
# Affiliation: Harvard University
# Contact: g.giecold@gmail.com; ggiecold@jimmy.harvard.edu


from __future__ import division

from operator import and_
from os import getcwd, makedirs, path
import re
from struct import calcsize, unpack
from sys import exit
from warnings import warn

import numpy as np
from rpy2.rinterface import NULL, TRUE
from rpy2.robjects import numpy2ri
from rpy2.robjects.packages import importr
from sklearn.manifold import TSNE

from . import Tree_classes


__all__ = ['Annotation', 'cytometry_preprocess', 'Cyto_data', 
           'FCS_handler', 'get_FCS_data', 'infer_pseudotime', 
           'PCR_preprocess', 'RNASeq_preprocess']


def infer_pseudotime(data, output_directory, tag = '', pcv_method = 'Rprincurve',
                     anchor_gene = None, markers = None):

    assert pcv_method in {'Rprincurve'} # taking into account the possibility of adding
                                        # in future versions other methods 
                                        # for principal curve analysis
    
    N_dim = 3
    model = TSNE(n_components = N_dim)
    TSNE_data = model.fit_transform(data)
    
    if pcv_method == 'Rprincurve':
        with open(path.join(output_directory, "{0}_TSNE_d{1}.tsv".format(tag, N_dim)),
                  'w') as f:
            f.write('\t'.join(['T{0}'.format(k) for k in xrange(1, N_dim + 1)]))
            f.write('\n')
            np.savetxt(f, TSNE_data, fmt = '%.6f', delimiter = '\t')
        
        numpy2ri.activate()
        princurve = importr('princurve')
        
        procedure = princurve.principal_curve
        fitpc = procedure(TSNE_data, NULL, 0.001, TRUE, 200, 2, 'lowess')
        curve_projections_matrix = np.array(fitpc.rx('s')[0])
        pseudotime_series = np.array(fitpc.rx('lambda')[0])
        
        with open(path.join(output_directory, "{0}_TSNE_d{1}_pcv.tsv".format(tag,
                  N_dim)), 'w') as f:
            np.savetxt(f, curve_projections_matrix, fmt = '%.6f', delimiter = '\t')
            
        with open(path.join(output_directory, "{0}_TSNE_d{1}_lambda.tsv".format(tag,
                  N_dim)), 'w') as f:
            np.savetxt(f, pseudotime_series, fmt = '%.6f', delimiter = '\t')
         
    else:
        print("ERROR: PySCUBA: Preprocessing: infer_pseudotime:\n"
              "your choice of method for principal curve analysis is not supported "
              "by the present version of PySCUBA.")
        exit(1)
        
    if anchor_gene:
        assert isinstance(anchor_gene, str)
        assert markers is not None
        
        N_cells_anchor = 1000
        
        gene_idx = np.where(markers == anchor_gene)[0]
        pseudotime_idx = np.argsort(pseudotime_series)
        
        anchor_gene_avg_beg = np.mean(data[pseudotime_idx[:N_cells_anchor], gene_idx])
        anchor_gene_avg_end = np.mean(data[pseudotime_idx[N_cells_anchor:], gene_idx])
        
        if anchor_gene_avg_end > anchor_gene_avg_beg:
            pseudotime_series = np.max(pseudotime_series) - pseudotime_series
        
    t_min = np.min(pseudotime_series)
    t_max = np.max(pseudotime_series)
    t_bins = 8
    
    cell_stages = t_bins * (pseudotime_series - t_min + 0.0001) / (t_max - t_min + 0.0002)
    cell_stages = np.ceil(cell_stages).astype(int).astype('str')
    
    return cell_stages
    

def parse_pairs(text):
    """Return (key, value) pairs from a string featuring particular delimiters.
       Modified from a corresponding function in the outdated 'fcm' 
       Python package by Jacob Frelinger.
    """
    
    delim = text[0]
    
    if delim == r'|':
        delim = '\|'
    elif delim == r'\a'[0]:
        delim = '\\\\'
        
    if delim != text[-1]:
        warn("WARNING: the text does not start and end with the same delimiter!")
        
    regex = re.compile('(?<=[^%s])%s(?!%s)' % (delim, delim, delim))
    tmp = text[1:-1].replace('$', '')
    tmp = regex.split(tmp)
    
    return dict(zip([x.lower() for x in tmp[::2]], tmp[1::2]))
    
    
class Annotation(object):
    """An annotation class instance stores meta-data from the recordings of a
       cytometry experiment.
       Modified from a corresponding class in the outdated 'fcm' 
       Python package by Jacob Frelinger.
    """

    def __init__(self, annotations = None):
        if annotations == None:
            annotations = {}

        self.__dict__['_mydict'] = annotations

    def __getattr__(self, name):
        if name in self._mydict.keys():
            self.__dict__[name] = self._mydict[name]
            
            return self._mydict[name]
        else:
            try:
                return self._mydict.__getattribute__(name)
            except:
                raise AttributeError("'{0}' has no attribue '{1}'".format(str(self.__class__), name))

    def __getstate__(self):
        return self._mydict

    def __setstate__(self, dict):
        self.__dict__['_mydict'] = dict
        
        for i in dict.keys():
            self.__dict__[i] = dict[i]

    def __setattr__(self, name, value):
        Annotation.__getattribute__(self, '_mydict')[name] = value
        self.__dict__[name] = value

    def __setitem__(self, name, value):
        self._mydict[name] = value
        self.__dict__[name] = self._mydict[name]

    def __getitem__(self, name):
        return self._mydict[name]

    def __repr__(self):
        return 'Annotation(' + self._mydict.__repr__() + ')'

    def __getstate__(self):
        return self.__dict__

    def __setstate(self, state):
        self.__dict__ = state

    def __getinitargs__(self):
        return (self._mydict,)

    def copy(self):
        return Annotation(self._mydict.copy())
        

class Cyto_data(object):
    """
    A Cyto_data object stores the data from a cytometry experiment.
    Modified from a corresponding class in the outdated 'fcm' Python 
    package by Jacob Frelinger.
    
    Members:
    --------
    Cyto_data.data_points : numpy.array
        The data points.
    
    Cyto_data.channels : list 
        Records which markers or scatters are in which columns.
    
    Cyto_data.scatters : list
        Keeps track of which indexes in Cyto_data.channels are scatters.
    """

    def __init__(self, name, data_points, channels, scatters = None, notes = None):
        """
        Parameters
        ----------
        name: name of the *.fcs file, barring any extension
        
        data_points: an array of data points
        
        channels: list
            Records which markers/scatters are in which columns.
            
        scatters: list 
            Which channels indexes denote scatters
        """
        
        self.name = name
        self.data_points = data_points
        self.tree = Tree_classes.Tree(data_points, channels)
        self.scatters = scatters
        
        self.markers = []
        if self.scatters is not None:
            for channel in range(len(channels)):
                if channel in self.scatters:
                    pass
                elif self.tree.root.channels[channel] in self.scatters:
                    pass
                else:
                    self.markers.append(channel)
                    
        if notes == None:
            notes = Annotation()
        self.notes = notes

    def __unicode__(self):
        return self.name

    def __repr__(self):
        return self.name

    def __getitem__(self, item):
        """Return the Cyto_data points.
        """
        
        if isinstance(item, tuple):
            item = list(item)
            if isinstance(item[1], str):
                item[1] = self.name_to_index(item[1])
            elif isinstance(item[1], tuple) or isinstance(item[1], list):
                item[1] = list(item[1])
                
                for i, j in enumerate(item[1]):
                    if isinstance(j, str):
                        print('{0} is string {1}'.format(i, j))
                        item[1][i] = self.name_to_index(j)
                        
            item = tuple(item)
                        
        return self.tree.view()[item]

    @property
    def channels(self):
        return self.current_node.channels
    
    def __getattr__(self, name):
            if name in dir(self.current_node.view()):
                return self.current_node.view().__getattribute__(name)
            else:
                raise AttributeError("'{0}' has no attribue '{1}'".format(str(self.__class__), name))

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, dict):
        for i in dict.keys():
            self.__dict__[i] = dict[i]

    def name_to_index(self, channels):
        """Return the channel indexes for the channels provided as arguments.
        """

        if isinstance(channels, str):
            return self.channels.index(channels)
            
        idx = []
        for i in channels:
            try:
                idx.append(self.channels.index(i))
            except ValueError:
                try:
                    for j in range(1, int(self.notes.text['par']) + 1):
                        if i == self.notes.text['p%dn' % j]:
                            idx.append(self.channels.index(self.notes.text['p%ds' % j]))
                except ValueError:
                    raise ValueError('{0} is not in list'.format(i))
                    
        if idx:
            return idx
        else:
            raise ValueError('The field {0} was not found'.format(str(channels)))

    def get_channel_by_name(self, channels):
        """Return the data associated with specific channel names.
        """
        
        return self.tree.view()[:, self.name_to_index(channels)]

    def get_markers(self):
        """Return the data associated with all the markers.
        """

        return self.view()[:, self.markers]

    def view(self):
        """Return the current view of the data.
        """
        
        return self.tree.view()

    def visit(self, name):
        """Switch the current view of the data.
        """
        
        self.tree.visit(name)

    @property
    def current_node(self):
        """Return the current node.
        """
        
        return self.tree.current

    def copy(self):
        """Return a copy of the Cyto_data object.
        """
        
        tname = self.name
        tree_data_points = self.tree.root.data
        tnotes = self.notes.copy()
        tchannels = self.channels[:]
        tscchannels = self.scatters[:]
        
        tmp = Cyto_data(tname, tree_data_points, tchannels, tscchannels, tnotes)
        
        from copy import deepcopy
        tmp.tree = deepcopy(self.tree)
        
        return tmp

    def gate(self, g, chan=None):
        """Return a gated region of the cytometry data.
        """
        
        return g.gate(self, chan)

    def subsample(self, s):
        """Return subsampled cytometry data.
        """
        
        return s.subsample(self)

    def get_current_node(self):
        """Return the current node.
        """
        
        return self.current_node

    def add_view(self, node):
        """Add a new node to the visualization tree.
        """
        
        self.tree.add_child(node.name, node)
        
        return self

    def summary(self):
        """Provide a summary of current view.
        """
        
        data_points = self.view()
        
        means = data_points.mean(0)
        stds = data_points.std(0)
        mins = data_points.min(0)
        maxs = data_points.max(0)
        medians = np.median(data_points, 0)
        
        dim = data_points.shape[1]
        
        summary = ''
        for i in range(dim):
            summary = summary + self.channels[i] + ":\n"
            summary = summary + "    max: " + str(maxs[i]) + "\n"
            summary = summary + "   mean: " + str(means[i]) + "\n"
            summary = summary + " median: " + str(medians[i]) + "\n"
            summary = summary + "    min: " + str(mins[i]) + "\n"
            summary = summary + "    std: " + str(stds[i]) + "\n"
            
        return summary

    def boundary_events(self):
        """Return a dictionary of the percentage of all events that are in the first 
           and in the last channel for each channel.
        """
        
        boundary_dict = {}
        for k, channel in enumerate(self.channels):
            col = self.view()[:, k]
            boundary_dict[channel] = \
                sum((col == min(col)) | (col == max(col))) / len(col)
                
        return boundary_dict
             
                
class NotSupportedFCSDataMode(Exception):
    """Exception raised for data modes in a *.fcs file that are not currently supported.
       Modified from a corresponding exception in the outdated 'fcm' 
       Python package by Jacob Frelinger.
    """

    def __init__(self, mode):
        self.mode = mode
        self.message = "Unfortunately, the FCS data stored as type {0} is not currently supported.".format(mode)
        self.args = (mode,)
        
        
def integer_format(b):
    """Return the binary format of an integer.
    """

    if b == 8:
        return 'B'
    elif b == 16:
        return 'H'
    elif b == 32:
        return 'I'
    else:
        print("Cannot handle integers of bit size {0}.".format(b))
        return None
        
        
def integer_bit_mask(b, ub):
    """Return the bit-mask of an integer and a bit-witdh.
    """

    if b == 8:
        return (0xFF >> (b - ub))
    elif b == 16:
        return (0xFFFF >> (b - ub))
    elif b == 32:
        return (0xFFFFFFFF >> (b - ub))
    else:
        print("Cannot handle integers of bit size {0}.".format(b))
        return None
        

def fluorescent_channel(name):
    """Check if a channel is a fluorescent channel.
    """
    
    name = name.lower()
    
    if name.startswith('cs'):
        return False
    elif name.startswith('fs'):
        return False
    elif name.startswith('ss'):
        return False
    elif name.startswith('ae'):
        return False
    elif name.startswith('cv'):
        return False
    elif name.startswith('time'):
        return False
    else:
        return True
        
                        
class FCS_handler(object):
    """Hold object to read and parse *.fcs files.
       Modified from a corresponding class in the outdated 'fcm' 
       Python package by Jacob Frelinger.
    """
    
    def __init__(self, file_path):
        self.file_path = file_path
        self.current_offset = 0

    def get_next_dataset(self, **kwargs):
        """Return the next cytometry dataset stored in a *.fcs file.
        """

        with open(self.file_path, 'rb') as self._f:
            header = self.parse_header(self.current_offset)

            text = self.parse_text(self.current_offset, header['text_start'], 
                                   header['text_stop'])
    
            try:
                analysis_beg = text['begin_analysis']
            except KeyError:
                analysis_beg = header['analysis_start']
                
            try:
                analysis_end = text['end_analysis']
            except KeyError:
                analysis_end = header['analysis_end']
                
            analysis = self.parse_analysis(self.current_offset, analysis_beg,
                                           analysis_end)
            
            try:
                data_beg = int(text['begin_data'])
            except KeyError:
                data_beg = header['data_start']
            try:
                data_end = int(text['end_data'])
            except KeyError:
                data_end = header['data_end']
    
            LMD = self.fix_LMD(self.current_offset, header['text_start'], 
                               header['text_stop'])
            data_end = data_end + LMD
            
            data = self.parse_data(self.current_offset, data_beg, data_end, text)
        
        channels = []
        scchannels = []
        scchannel_indexes = []
        base_chan_name = []
        
        for i in range(1, int(text['par']) + 1):    
            base_chan_name.append(text['p%dn' % i])
            
            try:
                if text['p%ds' % i] not in ['',' ']:
                    name = text['p%ds' % i]
                else:
                    name = text['p%dn' % i]
            except KeyError:
                name = text['p%dn' % i]
                
            channels.append(name)

            if not fluorescent_channel(name):
                scchannels.append(name)
                
                if name != 'Time':
                    scchannel_indexes.append(i - 1)
        
        _, name = path.split(self.file_path)
        name, _ = path.splitext(name)
        
        cyto_object = Cyto_data(name, data, channels, scchannels,
                                Annotation({'text': text,
                                            'header': header,
                                            'analysis': analysis,}))
                                                      
        return cyto_object


    def read_bytes(self, offset, start, stop):
        """Read bytes from start to stop, included.
        """

        self._f.seek(offset + start)

        return self._f.read(stop - start + 1)

    def parse_header(self, offset):
        """
        Parse the cytometry data in a *.fcs file at the specified offset 
        (accounting for the possibility of several data parts in the said file).
        """
        
        header = {}
        header['version'] = float(self.read_bytes(offset, 3, 5))
        header['text_start'] = int(self.read_bytes(offset, 10, 17))
        header['text_stop'] = int(self.read_bytes(offset, 18, 25))
        header['data_start'] = int(self.read_bytes(offset, 26, 33))
        header['data_end'] = int(self.read_bytes(offset, 34, 41))
        
        try:
            header['analysis_start'] = int(self.read_bytes(offset, 42, 49))
        except ValueError:
            header['analysis_start'] = -1
            
        try:
            header['analysis_end'] = int(self.read_bytes(offset, 50, 57))
        except ValueError:
            header['analysis_end'] = -1

        return header
        
    def parse_text(self, offset, start, stop):
        """Return the parsed text segment of a *.fcs file.
        """

        text = self.read_bytes(offset, start, stop)
        return parse_pairs(text)
        
    def parse_analysis(self, offset, start, stop):
        """Return the parsed analysis part of the *.fcs file under consideration.
        """

        if start == stop:
            return {}
        else:
            text = self.read_bytes(offset, start, stop)
            return parse_pairs(text)

    def fix_LMD(self, offset, start, stop):
        """Handle the LMD format (embedded FCS format) and the way it counts,
           which differs from other FCS formats.
        """
        
        text = self.read_bytes(offset, start, stop)
        
        if text[0] == text[-1]:
            return 0
        else:
            return -1

    def parse_data(self, offset, start, stop, text):
        """Return an array holding the data part of *.fcs file at hand.
        """

        dtype = text['datatype']
        mode = text['mode']
        tot = int(text['tot'])
        
        if mode == 'c' or mode == 'u':
            raise NotSupportedFCSDataMode(mode)

        if text['byteord'] == '1,2,3,4' or text['byteord'] == '1,2':
            order = '<'
        elif text['byteord'] == '4,3,2,1' or text['byteord'] == '2,1':
            order = '>'
        else:
            warn("WARNING: unsupported byte order {0}; using default @".format(text['byteord']))
            order = '@'

        bit_width = []
        data_range = []
        for i in range(1, int(text['par']) + 1):
            bit_width.append(int(text['p%db' % i]))
            data_range.append(int(text['p%dr' % i]))

        if dtype.lower() == 'i':
            data = self.parse_int_data(offset, start, stop, bit_width, data_range, 
                                       tot, order)
        elif dtype.lower() == 'f' or dtype.lower() == 'd':
            data = self.parse_float_data(offset, start, stop, dtype.lower(), tot, order)
        else:
            data = self.parse_ASCII_data(offset, start, stop, bit_width, dtype, 
                                         tot, order)
        
        return data

    def parse_int_data(self, offset, start, stop, bit_width, data_range, tot, order):
        """Parse *.fcs file and return data as an integer list.
        """

        if reduce(and_, [item in [8, 16, 32] for item in bit_width]):
            if len(set(bit_width)) == 1:
                num_items = (stop - start + 1) / calcsize(integer_format(bit_width[0]))
                tmp = unpack('%s%d%s' % (order, num_items, integer_format(bit_width[0])),
                             self.read_bytes(offset, start, stop))
            else:
                unused_bit_widths = map(int, map(np.log2, data_range))
                tmp = []
                
                current = start
                while current < stop:
                    for i, current_width in enumerate(bit_width):
                        bit_mask = integer_bit_mask(current_width, unused_bit_widths[i])
                        N_bytes = current_width / 8
                        bin_string = self.read_bytes(offset, current, current + N_bytes - 1)
                        current += N_bytes
                        val = bit_mask & unpack('%s%s' % (order, integer_format(current_width)), bin_string)[0]
                        tmp.append(val)
                        
        else:
            warn('WARNING: non-standard bit widths for the data part.')
            return None
            
        return np.array(tmp).reshape((tot, len(bit_width)))

    def parse_float_data(self, offset, start, stop, dtype, tot, order):
        """Parse a *.fcs file and return list of float data entries.
        """

        N_items = (stop - start + 1) / calcsize(dtype)
        tmp = unpack('%s%d%s' % (order, N_items, dtype), 
                     self.read_bytes(offset, start, stop))
                     
        return np.array(tmp).reshape((tot, len(tmp) / tot))

    def parse_ASCII_data(self, offset, start, stop, bit_width, dtype, tot, order):
        """Parse ASCII encoded data from a *.fcs file.
        """

        N_items = (stop - start + 1) / calcsize(dtype)
        tmp = unpack('%s%d%s' % (order, N_items, dtype), 
                     self.read_bytes(offset, start, stop))
        
        return np.array(tmp).reshape((tot, len(tmp) / tot))

    
def cytometry_preprocess(file_path, log_mode = False, pseudotime_mode = True,
                         pcv_method = 'Rprincurve', anchor_gene = None, 
                         exclude_marker_names = None):
                         
    data_tag, output_directory = create_output_directory(file_path)
                         
    cyto_object = get_FCS_data(file_path)
    
    marker_idx = np.array(cyto_objects.markers, dtype = str)
    marker_names = np.array(cyto_objects.channels[marker_idx], dtype = str)
    
    data = cyto_object.data_points
    data = data[:, marker_idx]
    
    cell_IDs = np.array(['cell_{0}'.format(i) for i in xrange(1, data.shape[0] + 1)], 
                        dtype = str)
    
    if exclude_marker_names:
        indices = np.zeros(0, dtype = int)
        for name in exclude_marker_names:
            indices = np.append(indices, np.where(marker_names == name)[0])
            
        data = np.delete(data, indices, axis = 1)
        marker_names = np.delete(marker_names, indices)
    
    cell_stages = infer_pseudotime(data, output_directory, data_tag, pcv_method,
                                   anchor_gene, marker_names)
                                   
    write_preprocessed_data(output_directory, cell_IDs, cell_stages, data, markers)
    
    return cell_IDs, data, marker_names, cell_stages.astype(float), data_tag, output_directory
            
      
def PCR_preprocess(file_path, log_mode = False, pseudotime_mode = False, 
                   pcv_method = 'Rprincurve', anchor_gene = None,
                   exclude_marker_names = None):
                   
    low_gene_fraction_max = 0.8
    
    data_tag, output_directory = create_output_directory(file_path)
    
    cell_IDs, cell_stages, data = get_PCR_or_RNASeq_data(file_path, pseudotime_mode)
    
    with open(file_path, 'r') as f:
        markers = np.loadtxt(f, dtype = str, delimiter = '\t', 
            skiprows = 1 if pseudotime_mode else 2, usecols = [0])
        markers.reshape(markers.size)
        
    if exclude_marker_names:
        indices = np.zeros(0, dtype = int)
        for name in exclude_marker_names:
            indices = np.append(indices, np.where(markers == name)[0])
            
        data = np.delete(data, indices, axis = 1)
        markers = np.delete(markers, indices)
        
    if pseudotime_mode:
        cell_stages = infer_pseudotime(data, output_directory, data_tag, pcv_method,
                                       anchor_gene, markers)
                                       
    condition = np.mean(data == 0, axis = 0) < low_gene_fraction_max
    data = np.compress(condition, data, 1)
    markers = np.compress(condition, markers)
    
    write_preprocessed_data(output_directory, cell_IDs, cell_stages, data, markers)
            
    return cell_IDs, data, markers, cell_stages.astype(float), data_tag, output_directory
 
                
def RNASeq_preprocess(file_path, log_mode = True, pseudotime_mode = False, 
                       pcv_method = 'Rprincurve', anchor_gene = None,
                       exclude_marker_names = None):

    assert isinstance(log_mode, bool)
    assert isinstance(pseudotime_mode, bool)
    
    # Threshold value for genes of low expression levels
    low_gene_threshold = 1
    # Maximum fraction of lowly-expressed cells allowed for each gene
    low_gene_fraction_max = 0.7
    # Number of highly variable genes selected
    N_selected_genes = 1000
    
    data_tag, output_directory = create_output_directory(file_path)
    
    cell_IDs, cell_stages, data = get_PCR_or_RNASeq_data(file_path, pseudotime_mode)
    
    with open(file_path, 'r') as f:
        markers = np.loadtxt(f, dtype = str, delimiter = '\t', 
            skiprows = 1 if pseudotime_mode else 2, usecols = [0])
        markers.reshape(markers.size)
        
    if exclude_marker_names:
        indices = np.zeros(0, dtype = int)
        for name in exclude_marker_names:
            indices = np.append(indices, np.where(markers == name)[0])
            
        data = np.delete(data, indices, axis = 1)
        markers = np.delete(markers, indices)
        
    if pseudotime_mode:
        cell_stages = infer_pseudotime(data, output_directory, data_tag, pcv_method,
                                       anchor_gene, markers)
        
    condition = np.mean(data < low_gene_threshold, axis = 0) < low_gene_fraction_max
    data = np.compress(condition, data, 1)
    markers = np.compress(condition, markers)
    
    Fano_factors = np.var(data, axis = 0) / np.mean(data, axis = 0).astype(float)
    idx = np.argsort(Fano_factors)[::-1][:N_selected_genes]
    data = data[:, idx]
    markers = markers[idx]
    
    if log_mode:
        np.log2(data + 1, data)
        
    write_preprocessed_data(output_directory, cell_IDs, cell_stages, data, markers)
        
    return cell_IDs, data, markers, cell_stages.astype(float), data_tag, output_directory
    

def create_output_directory(file_path):

    data_tag = path.basename(path.abspath(file_path)).split('.')[0]
    output_directory = path.join(getcwd(), 'SCUBA_analysis_of_{0}'.format(data_tag))
    try:
        makedirs(output_directory)
    except OSError:
        if not path.isdir(output_directory):
            raise
            
    return data_tag, output_directory
    
    
def get_FCS_data(file_path, **kwargs):
    """Return a data object from an *.fcs file"""

    cyto_object = FCS_handler(file_path)
    data = cyto_object.get_next_dataset(**kwargs)
    cyto_object._f.close()
    del cyto_object
    
    return data
    
    
def get_PCR_or_RNASeq_data(file_path, pseudotime_mode = False):

    with open(file_path, 'r') as f:
        cell_IDs = f.readline().rstrip('\n').split('\t')
        cell_IDs = np.array(cell_IDs[1:], dtype = str)
        
        if pseudotime_mode:
            cell_stages = np.empty(0, dtype = float)
        else:
            cell_stages = f.readline().rstrip('\n').split('\t')
            cell_stages = np.array(cell_stages[1:], dtype = str)
        
        data = np.loadtxt(f, dtype = float, delimiter = '\t', 
                          usecols = xrange(1, len(cell_IDs) + 1))     
        data = data.T
        
    return cell_IDs, cell_stages, data
    
    
def write_preprocessed_data(output_directory, cell_IDs, cell_stages, data, markers):

    processed_data_path = path.join(output_directory, 'processed_data.tsv')
    
    with open(processed_data_path, 'w') as f:
        f.write('\t'.join(cell_IDs))
        f.write('\n')
        f.write('\t'.join(cell_stages))
        f.write('\n')
        np.savetxt(f, data.T, fmt = '%.6f', delimiter = '\t')
        
    dataset = np.genfromtxt(processed_data_path, delimiter = '\t', dtype = str)
    dataset = np.insert(dataset, 0, np.append(['Cell ID', 'Stage'], 
        markers), axis = 1)
        
    with open(processed_data_path, 'w') as f:
        np.savetxt(f, dataset, fmt = '%s', delimiter = '\t')
        

    

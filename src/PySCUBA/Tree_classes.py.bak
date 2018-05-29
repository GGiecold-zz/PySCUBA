#!/usr/bin/env python


# PySCUBA/src/PySCUBA/Tree_classes.py


# Author: Gregory Giecold for the GC Yuan Lab
# Affiliation: Harvard University
# Contact: g.giecold@gmail.com; ggiecold@jimmy.harvard.edu


"""Classes for flow or mass cytometry data processing, modified from the outdated 'fcm' 
   Python package by Jacob Frelinger.
"""


import re

import numpy as np


__all__ = ['Tree']


class Node(object):
    """Base node object.
    """

    def __init__(self, name, parent, data):
        self.name = name
        self.parent = parent
        self.data = data
        self.prefix = 'n'

    def view(self):
        """Return the view of the data associated with this node.
        """

        return self.data

    def pprint(self, depth, size):
        tmp = "  " * depth + self.name
        
        if size:
            tmp = tmp + " " + str(self.view().shape[0])
    
        return tmp + "\n"

    def __getattr__(self, name):
        if name == 'channels':
            return self.parent.channels
        else:
            raise AttributeError("'{0}' has no attribue '{1}'".format(str(self.__class__), name))


class Root_node(Node):
    """Root node.
    """
    
    def __init__(self, name, data, channels):
        self.name = name
        self.parent = None
        self.data = data
        self.channels = channels
        self.prefix = 'root'


class Transform_node(Node):
    """Transformed data node.
    """

    def __init__(self, name, parent, data):
        self.name = name
        self.parent = parent
        self.data = data
        self.prefix = 't'

    def __getattr__(self, name):
        if name == 'channels':
            return self.parent.channels
        else:
            raise AttributeError("'{0}' has no attribue '{1}'".format(str(self.__class__), name))
        

class Subsample_node(Node):
    """Node of subsampled data.
    """

    def __init__(self, name, parent, param):
        self.name = name
        self.parent = parent
        self.param = param
        self.prefix = 's'
        
        if isinstance(param, tuple):
            self.channels = self.parent.channels[param[1]]

    def view(self):
        """Return the view of the data associated with this node.
        """
        
        return self.parent.view().__getitem__(self.param)


class Drop_channel_node(Node):
    """Node of data without some channels.
    """

    def __init__(self, name, parent, param, channels):
        self.name = name
        self.parent = parent
        self.param = param
        self.prefix = 'd'
        self.channels = channels

    def view(self):
        """Return the view of the data associated with this node.
        """
        
        return self.parent.view()[:, self.param]


class Gating_node(Node):
    """Node of gated data.
    """

    def __init__(self, name, parent, data):
        self.name = name
        self.parent = parent
        self.data = data
        self.prefix = 'g'

    def view(self):
        """Return the view of the data associated with this node.
        """
        
        if self.parent.view().shape[0] == 0:
            return np.array([]).reshape(self.parent.view().shape)
            
        return self.parent.view()[self.data]

    def __getattr__(self, name):
        if name == 'channels':
            return self.parent.channels
        else:
            raise AttributeError("'{0}' has no attribue '{1}'".format(str(self.__class__), name))


class Tree(object):
    """Tree of data for a Cyto_data object 
       (the latter is defined in the 'Preprocessing' module from the present package).
    """

    def __init__(self, data_points, channels):
        self.nodes = {}
        self.root = Root_node('root', data_points, channels)
        self.nodes['root'] = self.root
        self.current = self.root

    def parent(self):
        """Return the parent of a node.
        """
        
        return self.current.parent

    def children(self, node = None):
        """Return the children of a node."""
        
        if node == None:
            node = self.current
            
        return [i for i in self.nodes.values() if i.parent == node]

    def visit(self, name):
        """Visit a node in the tree."""

        if isinstance(name, str):
            self.current = self.nodes[name]
        elif isinstance(name, Node):
            self.current = name
        else:
            raise KeyError("No node named {0}.".format(name))

    def get(self, name = None):
        """Return the current node object."""
        
        if name is None:
            return self.current
        else:
            if name in self.nodes:
                return self.nodes[name]
            else:
                raise KeyError("No node named {0}.".format(name))

    def view(self):
        """Return a view of the current data.
        """
        
        return self.current.view()

    def add_child(self, name, node):
        """Add a node to the tree at the currently selected node.
        """
        
        if name == '':
            prefix = node.prefix
            pat = re.compile(prefix + "(\d+)")
            matches = [pat.search(i) for i in self.nodes]
            matches = [i for i in matches if i is not None]
            if len(matches):
                n = max([ int(i.group(1)) for i in matches])
                name = prefix + str(n + 1)
            else:
                name = prefix + '1'
                
        if name in self.nodes.keys():
            raise KeyError("Name {0} already in use in tree".format(name))
        else:
            node.name = name
            self.nodes[name] = node
            node.parent = self.current
            self.current = self.nodes[name]

    def rename_node(self, old_name, new_name):
        """Rename a node name;
           D(old,new) -> rename old to new.
        """
        
        if not self.nodes.has_key(old_name):
            raise KeyError("No node named {0}.".format(old_name))
        if self.nodes.has_key(new_name):
            raise KeyError("There is already a node named {0}".format(new_name))
        else:
            self.nodes[new_name] = self.nodes[old_name]
            self.nodes[new_name].name = new_name
            del self.nodes[old_name]

    def pprint(self, size = False):
        return self._rpprint(self.root, 0, size)

    def _rpprint(self, n, d, size = False):
        tmp = n.pprint(d, size)
        for i in self.children(n):
            tmp += self._rpprint(i, d + 1, size)
            
        return tmp




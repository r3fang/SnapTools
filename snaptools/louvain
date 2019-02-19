# -*- coding: utf-8 -*-
""" 

The MIT License

Copyright (c) 2018 Rongxin Fang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

"""

import networkx as nx
import community
from argparse import ArgumentParser, RawTextHelpFormatter
import os

def louvain(
    edge_file,
    output_file,
    resolution
    ):
    
    ################################################################################################      
    if not os.path.exists(edge_file):
        print('Error: ' + edge_file + ' does not exist!');
        sys.exit(1);
	
    ################################################################################################      
    G = nx.Graph();
    with open(edge_file) as fin:
        for line in fin:
            elems = line.split();
            G.add_edge(elems[0],elems[1], weight=float(elems[2]))
    
    partition = community.best_partition(G, resolution = resolution)        
    
    with open(output_file, "w") as fout:
        for node in partition.keys():
            fout.write(str(node) + " " + str(partition[node] + 1) + "\n") 


'''
The following program computes the topological properties of a cloud of points using a coordinates file as input
We must have in the same folder where the program is executed a 'points.txt' file with numbers stored in the following format:
       1.23 4.56
       6.78 9.01
       12.43 -12.02
They can be either integer or float

Additionally, we must have the 'TDA.py' file with the Topological analysis functions created at:
    https://github.com/multinetlab-amsterdam/network_TDA_tutorial

As an output, it creates a 'Results' folder with 4 new files:
    - Result_SC.txt: with the abstract information of the resulting simplicial complex
    - Geometrical_SC.png: geometrical representation of the simplicial complex, highlighting the relevant simplices
    - Persistence_barcode.png: plot of the persistence homology properties
'''


# Import required Python packages
import numpy as np # version 1.18.5
import networkx as nx # version 2.4
import community # version 0.13 (python-louvain)
import gudhi # version 3.3.0
import itertools

import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import scipy.io # version 1.4.1
from sklearn import preprocessing # version 0.23.1
import itertools

import seaborn as sns
import pandas as pd
import glob
from pathlib import Path
import time

get_ipython().run_line_magic('load_ext', 'watermark')

import Simplicial_Complexes as sc

import os
if not os.path.exists('Results'):
    os.makedirs('Results')

# Create Graph class
class Graph(object):

    # Initialize the matrix
    def __init__(self, size):
        self.adjMatrix = []
        for i in range(size):
            self.adjMatrix.append([0 for i in range(size)])
        self.size = size

    # Add edges
    def add_edge(self, v1, v2):
        if v1 == v2:
            print("Same vertex %d and %d" % (v1, v2))
        self.adjMatrix[v1][v2] = 1
        self.adjMatrix[v2][v1] = 1

    # Remove edges
    def remove_edge(self, v1, v2):
        if self.adjMatrix[v1][v2] == 0:
            print("No edge between %d and %d" % (v1, v2))
            return
        self.adjMatrix[v1][v2] = 0
        self.adjMatrix[v2][v1] = 0

    def __len__(self):
        return self.size

    # Print the matrix
    def print_matrix(self):
        print(np.array(self.adjMatrix))
        
# Load the set of points (input raw data)
data = list()
for line in open('points.txt', 'r'): 
    row = line.split()
    data.append((float(row[0]),float(row[1])))    

# Parameters
parameters_list = list()
for line in open('parameters.txt', 'r'):
    row = line.split()
    parameters_list.append(float(row[0]))
    
field_extension = parameters_list[0] #Maximum entension of a synaptic field to create an edge
min_d = int(parameters_list[1])
max_d = int(parameters_list[2])

maxsimplex = max_d #We need to specify a maximum dimension of a simplex to make a computationally efficient program
n_vertices = len(data)
position_dict = {i:data[i] for i in range(n_vertices)} # To specify the position of our nodes in networkx package

points = np.asarray(data)
g = Graph(n_vertices)

# Create simplicial complex
start_time = time.time()
cech_complex = sc.EnumSimplices(points, field_extension, max_d)
print('Number of points: ', n_vertices)
print('Maximum dimension of simplices: ', max_d)
print("--- %s seconds ---" % (time.time() - start_time))
#cech_complex = sc.EnumMaxSimplices_Efficient(points,field_extension,list(),list(np.arange(n_vertices)))
st = gudhi.SimplexTree()
for simplex in cech_complex:
    st.insert(simplex, float(len(simplex)))

# Note that, if there are n vertices, then the maximum number of edges is n*(n-1)
result_str = 'Simplicial complex is of dimension ' + repr(st.dimension()) + ' - ' +     repr(st.num_simplices()) + ' simplices - ' +     repr(n_vertices) + ' vertices.'

filtration = st.get_filtration()

with open('./Results/Result_SC.txt','w+') as file:
    file.write(result_str)
    file.write('\nAbstract simplicial complex: (Simplex -> Filtration value)\n')
    for filtered_value in filtration:
        file.write("%s -> %.2f\n" % tuple(filtered_value))
        simplex = filtered_value[0]
        if(len(simplex)==2):
            g.add_edge(simplex[0],simplex[1])
    
cmap = cm.hot
norm = Normalize(vmin=0, vmax=max_d)
m = cm.ScalarMappable(norm=norm, cmap=cmap)

# Visualize matrix as graph, highlighting the relevant simplices
A = np.array(g.adjMatrix)
G = nx.from_numpy_matrix(A)  

plt.figure(figsize=(20,16))
pos = position_dict  # positions for all nodes

options = {"node_color": "black", "edgecolors": "black", "alpha": 0.9}

nx.draw_networkx_nodes(G, pos, **options)
nx.draw_networkx_edges(G, pos, edge_color='k') # Background network with all edges




filtration = st.get_filtration()
for filtered_value in filtration:
    simplex = filtered_value[0]
    if(sc.IsMaxCell(points, np.array(simplex).astype(int), field_extension)):
        dim = len(simplex)
        if dim >= min_d:
            edges = list(itertools.combinations(simplex, 2)) #Return 2 length subsequences of elements from the input simplex.
            nx.draw_networkx_edges(
                G,
                pos,
                edgelist=edges,
                width=1.5,
                edge_color = 'r'
            )

labels = {i:i for i in range(n_vertices)}

nx.draw_networkx_labels(G, pos, labels, font_color="whitesmoke")

plt.title('Geometrical simplicial complex')
#plt.grid()
#plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap.reversed()))
plt.draw()
plt.savefig('./Results/Geometrical_SC.png')

plt.figure()
diag = st.persistence(min_persistence=-1)
gudhi.plot_persistence_barcode(diag, legend=True, max_intervals=0)
plt.savefig('./Results/Persistence_barcode.png')
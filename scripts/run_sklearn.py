import numpy as np
from sklearn.neighbors import NearestNeighbors
import sys

# Usage:
# - From command line: 
#     python run_sklearn.py path_to_input_data_file k d N algorithm
#     example: python run_sklearn.py ../data/frey.dat 12 560 1965 kd_tree sklearn_output.txt
# - In interpreter: 
#     knn_graph = run_sklearn.get_knn_graph(path_to_input, k, d, N, algorithm)
#     example: knn_graph = run_sklearn.get_knn_graph('../data/frey.dat', 12, 560, 1965, 'kd_tree')

def get_knn_graph(data_file, k, d, N, alg):
  a = np.fromfile(data_file, dtype=float).reshape((N,d))
  k_plus_1 = k+1
  nbrs = NearestNeighbors(n_neighbors=(k_plus_1), algorithm=alg).fit(a)
  knn_graph = nbrs.kneighbors_graph(a).toarray()
  return knn_graph

if __name__ == '__main__':
  knn_graph = get_knn_graph(sys.argv[1:], int(argv[1]), int(argv[2]), int(argv[3]), argv[4])
  knn_graph.tofile(sys.argv[6], "\t")
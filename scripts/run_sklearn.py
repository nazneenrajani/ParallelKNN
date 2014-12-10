import numpy as np
from sklearn.neighbors import NearestNeighbors
import time
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
  t_start = time.time()
  nbrs = NearestNeighbors(n_neighbors=(k_plus_1), algorithm=alg).fit(a)
  t_tree = time.time()
  knn_graph = nbrs.kneighbors_graph(a).toarray()
  t_graph = time.time() - t_tree
  t = time.time() - t_start
  print 'time to fit model = ' + str(t_tree-t_start) + " seconds"
  print 'time to make knn graph = ' + str(t_graph) + " seconds"
  print 'overall time = ' + str(t) + " seconds"
  return knn_graph

if __name__ == '__main__':
  knn_graph = get_knn_graph(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
  # knn_graph.tofile(sys.argv[6], "\t")
import numpy as np
from sklearn.neighbors import NearestNeighbors

# from command line: python run_sklearn.py chen_knn/example.dat 12 560 1965
# in interpreter: 

def get_knn_graph(argv):
  data_file = argv[0]
  k = argv[1]
  d = argv[2]
  N = argv[3]
  a = np.fromfile(data_file, dtype=float).reshape((N,d))
  k_plus_1 = k+1
  nbrs = NearestNeighbors(n_neighbors=(k_plus_1), algorithm='kd_tree').fit(a)
  knn_graph = nbrs.kneighbors_graph(a).toarray()
  return knn_graph

if __name__ == 'main':
  get_knn_graph(sys.argv[1:])
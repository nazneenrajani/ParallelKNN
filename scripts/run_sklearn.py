import numpy as np
from sklearn.neighbors import NearestNeighbors
import sys

# from command line: python run_sklearn.py chen_knn/example.dat 12 560 1965 kd_tree sklearn_output.txt
# in interpreter: knn_graph = run_sklearn.get_knn_graph(['chen_knn/example.dat', 12, 560, 1965, 'kd_tree'])

def get_knn_graph(argv):
  data_file = argv[0]
  k = int(argv[1])
  d = int(argv[2])
  N = int(argv[3])
  alg = argv[4]
  a = np.fromfile(data_file, dtype=float).reshape((N,d))
  k_plus_1 = k+1
  nbrs = NearestNeighbors(n_neighbors=(k_plus_1), algorithm=alg).fit(a)
  knn_graph = nbrs.kneighbors_graph(a).toarray()
  return knn_graph

if __name__ == '__main__':
  knn_graph = get_knn_graph(sys.argv[1:])
  knn_graph.tofile(sys.argv[6], "\t")
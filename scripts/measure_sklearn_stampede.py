import sys
import run_sklearn

k = 8
d = 784
datasets = [['../../data/mnist-test.dat', 10000], ['../../data/mnist-train.dat', 60000]]
for dataset in datasets:
  print dataset[0], ' -----------------------'
  algs = ['kd_tree', 'ball_tree', 'brute']
  for alg in algs:
    print alg, ": +++++++++++++"
    run_sklearn.get_knn_graph(dataset[0], k, d, dataset[1], alg)
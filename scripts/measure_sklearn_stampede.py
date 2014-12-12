import sys
import run_sklearn

poker = ["/scratch/02234/kmcardle/data/poker.t", 'libsvm', 10, 1000000, 8]
rna = ["/scratch/02234/kmcardle/data/cod-rna.t", 'libsvm', 8, 271617, 8]
cadata = ["/scratch/02234/kmcardle/data/cadata", 'libsvm', 8, 20640, 8]
covtype = ["/scratch/02234/kmcardle/data/covtype.libsvm.binary", 'libsvm', 54, 581012, 8]
year = ["/scratch/02234/kmcardle/data/YearPredictionMSD", 'libsvm', 90, 463715, 8]
aloi = ["/scratch/02234/kmcardle/data/aloi", 'libsvm', 128, 108000, 8]

# datasets = [poker, rna, cadata, covtype, year, aloi]
datasets = [rna, covtype, poker, year, aloi]

for dataset in datasets:
  print dataset[0], ' -----------------------'
  algs = ['kd_tree', 'ball_tree', 'brute']
  for alg in algs:
    print alg, ": +++++++++++++"
    # get_knn_graph(data_file, data_format, k, d, N, alg):
    run_sklearn.get_knn_graph(dataset[0], dataset[1], dataset[4], dataset[2], dataset[3], alg)
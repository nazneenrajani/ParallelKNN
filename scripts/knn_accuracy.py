import numpy as np
import sys
import run_sklearn

# from command line: python knn_accuracy.py chen_knn/example.dat chen_knn/example.knn chen 12 560 1965
# in interpreter: knn_accuracy.get_accuracy(['chen_knn/example.dat', 'chen_knn/example.knn', 'chen', 12, 560, 1965)

def get_accuracy(argv):
  input_data = argv[0]
  knn_output = argv[1]
  alg = argv[2]
  k = int(argv[3])
  d = int(argv[4])
  N = int(argv[5])

  # Get ground truth
  ground_truth = run_sklearn.get_knn_graph([input_data, k, d, N, 'brute'])
  
  if alg == "chen":
    # Get knn graph file in numpy array
    knn_graph = np.zeros((N,N), dtype=np.int)
    with open(knn_output, 'r') as f:
      for line in f:
        line.strip()
        line = line.split()
        j = int(line[0]) - 1
        i = int(line[1]) - 1
        knn_graph[i][j] = 1

  elif alg == "sklearn":
    knn_graph = knn_output

  if (knn_graph == ground_truth).all():
    print "100% accuracy"
  else:
    print "not 100% accurate"

if __name__ == '__main__':
  get_accuracy(sys.argv[1:])
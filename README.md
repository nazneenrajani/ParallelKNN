ParallelKNN
===========

Parallel implementation of KNN graph using KD Trees / Ball Trees

This repository is the final project for Scalable Machine Learning. 4 approaches to solving the parallel KNN problem are provided.

1. KDTrees with OpenMP : Uses OpenMP to parallelize KNN construction using KDTrees
2. KDTrees with Galois : Uses Galois to construct and parallelize KNN construction using KDTrees
3. Ball Trees with OpenMP : Uses OpenMP to parallelize KNN construction using Ball Trees
4. KDTrees with Galois : Uses Galois to construct and parallelize KNN construction using KDTrees

All experiments are run on Stampede (part of TACC). The baseline is the brute force construction of KNN graph.

Please refer to  the report for indepth details of each of the approaches. The report and data sets can be found at http://www.cs.utexas.edu/~nrajani/index.html#project

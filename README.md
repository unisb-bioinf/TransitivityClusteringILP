# TransitivityClusteringILP
Clustering algorithm based on the clique partitioning problem.

C++ implementation of the Integer Linear Programming approach of Grötschel and Wakabayashi [1].
In out implementation, we use slightly different constaints that allow us to calculate clusters itaritively.

# Requirements 

C++ Compiler (C++14 ready)
Boost (Version 1.6 or higher)
Eigen (Version 4 or higher)
CPLEX (Version 1262)
GTEST

# Installation

mkdir build
cd build
cmake .. -DCPLEX_SRC_DIR=$PATH
make

# Usage

-m [ --matrix ] arg                                     Name of a text file containing a tab-separated matrix.
                                                        Example:
                                                        Sample1 Sample2 Sample3
                                                        Feature1 1.0  2.0 3.0
                                                        Feature2 2.4  0.0 1.4
                                                        ...

-r [ --no-row-names ]                                   Indicator that the file does not contain row names. (optional)
-c [ --no-col-names ]                                   Indicator that the file does not contain column names. (optional)
-t [ --cplex-threads ] arg (=32)                        Max number of threads that should be used by CPLEX. (optional)
-l [ --cplex-time-limit ] arg (=10)                     Time limit (in min) for each CPLEX iteration. (optional)
-d [ --dissimilarity ]                                  Indicator that the matrix already is a dissimilarity matrix.
-s [ --similarity-measure ] arg (=pearson-correlation)  Method used to compute pairwise dissimilarities.
                                                         Measures:
                                                         pearson-correlation
                                                         spearman-correlation
                                                         distance-correlation
                                                         euclidean-distance
                                                         euclidean-distance-for-gradients
                                                         shifted-euclidean-distance-for-points
                                                         angle-distance
                                                         dynamic-time-warping

-i [ --run-itratively ]                                  Run algorithm iteratively util threshold, maximum number of iteration, or time limit is reached.
-f [ --fix-previous-solution ]                           Fix the solution of the previous iteration. (This makes the algorithm equivalent to hierachical clustering with complete linkage)
-t [ --similarity-threshold ] arg (=0.94999999999999996) Threshold for minimum dissimilarity score that should be used.
-x [ --max-number-of-iterations ] arg (=1000)            Max number of iterations.
-k [ --skip-k-iterations ] arg (=0)                      Number of iterations that should be performed before first model should be solved.
-p [ --output-prefix ] arg                               Prefix for all intermediate output files.
-e [ --output-edges ] arg                                Output file for sorted edge.
-o [ --output ] arg                                      Name of the final output file.

# References

[1] Grötschel, M., Wakabayashi, Y. A cutting plane algorithm for a clustering problem. Mathematical Programming 45, 59–96 (1989). https://doi.org/10.1007/BF01589097

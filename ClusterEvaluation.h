/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2019 Tim Kehl tkehl@bioinf.uni-sb.de>
 *
 */
#ifndef TRANSITIVITY_CLUSTERING_ILP_CLUSTER_EVALUATION_H
#define TRANSITIVITY_CLUSTERING_ILP_CLUSTER_EVALUATION_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>
#include <iostream>
#include <map>
#include <set>

#include <boost/math/special_functions/atanh.hpp>

#include "macros.h"

#include "DenseMatrix.h"

namespace TransitivityClusteringILP
{

/**
 * A collection of mathematical operations.
 */
namespace ClusterEvaluation
{

/**
 * This method extracts the clusters from an ILP result.
 *
 * @param matrix Binary adjacency matrix
 *
 * @return clusters
 */
template<typename Matrix>
std::vector<std::vector<size_t>> extract_clusters(const Matrix& matrix) {
	std::vector<std::vector<size_t>> clusters;
	std::set<size_t> skip;
	std::set<size_t>::iterator it;
	for(size_t i=0; i<matrix.rows(); ++i) {
		it = skip.find(i);
		if(it != skip.end()) {
			continue;
		}
		std::vector<size_t> cluster;
		for(size_t j=i; j<matrix.cols(); ++j) {
			if(matrix(i,j) == 1.0) {
				cluster.emplace_back(j);
				skip.emplace(j);
				skip.emplace(i);
			}
		}
		clusters.emplace_back(cluster);
	}
	return std::move(clusters);
}

/**
 * This method calculates centroids of the given clusters.
 *
 * @param expression Expression matrix
 * @param clusters Clustering
 *
 * @return Centroids of clusters
 */
/*template<typename Matrix>
std::vector<std::vector<double>> centroids(const Matrix& expression, const std::vector<std::vector<size_t>>& clusters) {
	std::vector<std::vector<double>> centroids(expression.cols(), std::vector<double>(expression.cols(), 0.0));
	double n = (double) clusters.size();
	for (size_t i=0; i<clusters.size(); ++i) {
		for(size_t j=0; j<clusters[i].size(); ++j) {
			for (size_t k=0; k<expression.cols(); ++k) {
				centroids[i][k] += expression(clusters[i][j], k);
			}
		}	
	}
	for (size_t i=0; i<centroids.size(); ++i) {
		centroids[i] /= n;
	}
	return std::move(centroids);
}
*/
/**
 * This method calculates the Davies-Bouldin index of a given clustering.
 *
 * Since algorithms that produce clusters with low intra-cluster distances (high intra-cluster similarity) 
 * and high inter-cluster distances (low inter-cluster similarity) will have a low Davies–Bouldin index, 
 * the clustering algorithm that produces a collection of clusters with the smallest Davies–Bouldin index 
 * is considered the best algorithm based on this criterion.
 *
 * @param meas Similarity measure
 * @param expression Expression matrix
 * @param clusters Clusters
 *
 * @return Davies-Bouldin index
 */
/*template<typename SimilarityMeasure, typename Matrix>
double davies_bouldin_index(SimilarityMeasure meas, const Matrix& expression, const std::vector<std::vector<size_t>>& clusters) {
	double db = 0.0;
	double n = (double)dist.rows();
	std::vector<std::vector<double>> centroids = centroids(expression, clusters);
	// Calculate intra-cluster differences
	std::vector<double> sigmas(clusters.size(), 0.0);
	for(size_t i=0; i<clusters.size(); ++i) {
		for(size_t j=0; j<clusters[i].size(); ++j) {
			RowMajorMatrixIterator<Matrix> rit(&expression, clusters[i][j]);
			std::vector<double> gene(rit->begin(), rit->end());
			// Since we use a similarity measure instead of a distance measure,
			// we use d = 1 - similatiry
		        double d = 1.0 - meas.compute_similarity(
		        	gene.begin(),
		                gene.end(),
		                centroids[i].begin(),
		                centroids[i].end()
		        );
			sigmas[i] += d;
		}
		sigmas[i] /= (double) clusters[i].size();
	}

	// Calculate differences between clusters 
	for(size_t i=0; i<centroids.size(); ++i) {
		double max = 0.0;
		for(size_t j=0; j<centroids.size(); ++j) {
			if(j == i) continue;
			// Since we use a similarity measure instead of a distance measure,
			// we use d = 1 - similatiry
			double d = 1.0 - meas.compute_similarity(
		        	centroids[i].begin(),
		                centroids[i].end(),
		                centroids[j].begin(),
		                centroids[j].end()
		        );
			double dbi = (sigmas[i] + sigmas[j]) / d;
			max = std::max(max, dbi);
		}
		db += max;
	}
	
	return db / n;
}*/

/**
 * This method calculates the Dunn index of a given clustering.
 *
 * Since internal criterion seek clusters with high intra-cluster similarity and low inter-cluster similarity, 
 * algorithms that produce clusters with high Dunn index are more desirable.
 *
 * @param meas Similarity measure
 * @param expression Expression matrix
 * @param clusters Clusters
 *
 * @return Davies-Bouldin index
 */
/*template<typename SimilarityMeasure, typename Matrix>
double dunn_index(SimilarityMeasure meas, const Matrix& expression, const std::vector<std::vector<size_t>>& clusters) {
	double d = 0.0;
	std::vector<std::vector<double>> centroids = centroids(expression, clusters);

	// Calculate minimum distance between clusters
	double dij_min = 0.0;
	for(size_t i=0; i<clusters.size(); ++i) {
		for(size_t j=i+1; j<clusters.size(); ++j) {
			if(j == i) continue;
			// Since we use a similarity measure instead of a distance measure,
			// we use d = 1 - similatiry
			double d = 1.0 - meas.compute_similarity(
		        	centroids[i].begin(),
		                centroids[i].end(),
		                centroids[j].begin(),
		                centroids[j].end()
		        );
			dij_min = std::min(d, dij_min);
		}
	}	

	// Calculate maximum intra-cluster distances
	double dk_max = 0.0;
	// Calculate intra-cluster differences
	for(size_t i=0; i<clusters.size(); ++i) {
		double di = 0.0;
		for(size_t j=0; j<clusters[i].size(); ++j) {
			RowMajorMatrixIterator<Matrix> rit(&expression, clusters[i][j]);
			std::vector<double> gene(rit->begin(), rit->end());
			// Since we use a similarity measure instead of a distance measure,
			// we use d = 1 - similatiry
		        double d = 1.0 - meas.compute_similarity(
		        	gene.begin(),
		                gene.end(),
		                centroids[i].begin(),
		                centroids[i].end()
		        );
			di += d;
		}
		di /= (double) clusters[i].size();
		dk_max = std::max(dk_max, di);
	}

	d = dij_min / dk_max;

	return d;
}*/

}
}

#endif // TRANSITIVITY_CLUSTERING_ILP_CLUSTER_EVALUATION_H

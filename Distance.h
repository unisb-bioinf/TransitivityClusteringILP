/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2019 Tim Kehl tkehl@bioinf.uni-sb.de>
 *
 */
#ifndef TRANSITIVITY_CLUSTERING_ILP_STATISTIC_H
#define TRANSITIVITY_CLUSTERING_ILP_STATISTIC_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>
#include <iostream>
#include <map>

#include <boost/math/special_functions/atanh.hpp>

#include "macros.h"

#include "DenseMatrix.h"

namespace TransitivityClusteringILP
{

/**
 * A collection of mathematical operations.
 */
namespace Distance
{

/**
 * This methods implements the squared distance.
 *
 * @param first_begin  InputIterator corresponding to the start of the first group.
 * @param first_end    InputIterator corresponding to the end of the first group.
 * @param second_begin InputIterator corresponding to the start of the second group.
 * @param second_end   InputIterator corresponding to the end of the second group.
 *
 * @return Distance
 */
template <typename value_type, typename InputIterator>
value_type
mean_difference(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator)
{
	value_type n = (value_type)std::distance(first_begin, first_end);
	value_type diff = 0;
	for(; first_begin != first_end; ++first_begin, ++second_begin) {
		diff += (*first_begin - *second_begin);
	}
	return diff / n;
}

/**
 * This methods implements the squared distance.
 *
 * @param first_begin  InputIterator corresponding to the start of the first group.
 * @param first_end    InputIterator corresponding to the end of the first group.
 * @param second_begin InputIterator corresponding to the start of the second group.
 * @param second_end   InputIterator corresponding to the end of the second group.
 *
 * @return Distance
 */
template <typename value_type, typename InputIterator>
value_type
squared_distance(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator)
{
	value_type n = (value_type)std::distance(first_begin, first_end);
	value_type dist = 0;
	for(; first_begin != first_end; ++first_begin, ++second_begin) {
		value_type diff = (*first_begin - *second_begin);
		dist += diff * diff;
	}
	return dist / n;
}

/**
 * This methods implements the euclidean distance.
 *
 * @param first_begin  InputIterator corresponding to the start of the first group.
 * @param first_end    InputIterator corresponding to the end of the first group.
 * @param second_begin InputIterator corresponding to the start of the second group.
 * @param second_end   InputIterator corresponding to the end of the second group.
 *
 * @return Distance
 */
template <typename value_type, typename InputIterator>
value_type
euclidean_distance(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator second_end)
{
	return std::sqrt(squared_distance<value_type, InputIterator>(first_begin, first_end, second_begin, second_end));
}

/**
 * This methods implements the shifted euclidean distance.
 *
 * @param first_begin  InputIterator corresponding to the start of the first group.
 * @param first_end    InputIterator corresponding to the end of the first group.
 * @param second_begin InputIterator corresponding to the start of the second group.
 * @param second_end   InputIterator corresponding to the end of the second group.
 *
 * @return Distance
 */
template <typename value_type, typename InputIterator>
value_type
shifted_euclidean_distance_for_points(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator)
{
	value_type n = (value_type)std::distance(first_begin, first_end);
	value_type s = 0.0;
	// Calculate optimal shift
	InputIterator first = first_begin;
	InputIterator second = second_begin;
	for(; first != first_end; ++first, ++second) {
		s += (*first - *second);
	}
	s = s / n;

	value_type dist = 0.0;
	//Calculate shifted distance
	for(; first_begin != first_end; ++first_begin, ++second_begin) {
		value_type diff = (*first_begin - s - *second_begin);
		dist += diff * diff;
	}
	return std::sqrt(dist / n);
}

/**
 * This methods implements the shifted euclidean distance.
 *
 * @param first_begin  InputIterator corresponding to the start of the first group.
 * @param first_end    InputIterator corresponding to the end of the first group.
 * @param second_begin InputIterator corresponding to the start of the second group.
 * @param second_end   InputIterator corresponding to the end of the second group.
 *
 * @return Distance
 */
template <typename value_type, typename InputIterator>
value_type
shifted_euclidean_distance_for_gradients(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator second_end)
{
	std::vector<value_type> fc_1;
	fc_1.reserve(std::distance(first_begin, first_end));
	for(; first_begin != (first_end-1); ++first_begin) {
		fc_1.emplace_back(*(first_begin+1) - *first_begin);
	}
	std::vector<value_type> fc_2;
	fc_2.reserve(std::distance(second_begin, second_end));
	for(; second_begin != (second_end-1); ++second_begin) {
		fc_2.emplace_back(*(second_begin+1) - *second_begin);
	}
	return shifted_euclidean_distance_for_points<value_type>(fc_1.begin(), fc_1.end(), fc_2.begin(), fc_2.end());
}

/**
 * This methods implements the shifted euclidean distance.
 *
 * @param first_begin  InputIterator corresponding to the start of the first group.
 * @param first_end    InputIterator corresponding to the end of the first group.
 * @param second_begin InputIterator corresponding to the start of the second group.
 * @param second_end   InputIterator corresponding to the end of the second group.
 *
 * @return Distance
 */
template <typename value_type, typename InputIterator>
value_type
shifted_euclidean_distance_for_gradients_and_points(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator second_end)
{

	return 	shifted_euclidean_distance_for_points<value_type>(first_begin, first_end, second_begin, second_end) +
		shifted_euclidean_distance_for_gradients<value_type>(first_begin, first_end, second_begin, second_end);
}

template <typename value_type, typename InputIterator>
GeneTrail::DenseMatrix pairwise_distance(InputIterator begin, InputIterator end) {
	std::vector<value_type> x(begin, end);
	GeneTrail::DenseMatrix dist(x.size(), x.size());
	for(size_t i=0; i<x.size(); ++i){
                dist.set(i, i, 0.0);
                for(size_t j=i+1; j<x.size(); ++j){
			value_type diff = x[i] - x[j];
			dist(i, j) = std::sqrt(diff*diff);
			dist(j, i) = dist(i, j);
		}
	}
	return std::move(dist);
}

GeneTrail::DenseMatrix normalize_distances(GeneTrail::DenseMatrix& A) {
	double n = (double)A.rows();
	double mean = 0.0;
	std::vector<double> row_means(A.rows(), 0.0);	
	std::vector<double> col_means(A.cols(), 0.0);
	for(size_t i=0; i<A.rows(); ++i){
		for(size_t j=0; j<A.cols(); ++j){
			mean += A(i,j);
			row_means[i] += A(i,j);
			col_means[j] += A(i,j);
		}
	}
	mean /= n*n;
	for(size_t i=0; i<n; ++i){
		row_means[i] /= n;
		col_means[i] /= n;
	}
	for(size_t i=0; i<A.rows(); ++i){
		for(size_t j=0; j<A.cols(); ++j){
			A(i,j) = A(i,j) - row_means[i] - col_means[j] + mean;
		}
	}
	return A;
}

/**
 * This methods implements the distance covariance.
 *
 * @return Distance
 */
template <typename value_type>
value_type distance_covariance(GeneTrail::DenseMatrix& A, GeneTrail::DenseMatrix& B)
{
	double n = (double)A.rows();
	value_type cov = 0.0;
	for(size_t i=0; i<A.rows(); ++i){
		for(size_t j=0; j<A.cols(); ++j){
			cov += A(i,j)*B(i,j);
		}
	}
	cov /= n*n;
	return(std::sqrt(cov));
}

/**
 * This methods implements the distance correlation.
 *
 * @param first_begin  InputIterator corresponding to the start of the first group.
 * @param first_end    InputIterator corresponding to the end of the first group.
 * @param second_begin InputIterator corresponding to the start of the second group.
 * @param second_end   InputIterator corresponding to the end of the second group.
 *
 * @return Distance
 */
template <typename value_type, typename InputIterator>
value_type distance_covariance(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator second_end)
{
	GeneTrail::DenseMatrix A = pairwise_distance<value_type>(first_begin, first_end);
	A = normalize_distances(A);
	GeneTrail::DenseMatrix B = pairwise_distance<value_type>(second_begin, second_end);
	B = normalize_distances(B);
	return distance_covariance<value_type>(A, B);
}

/**
 * This methods implements the distance correlation.
 *
 * @param first_begin  InputIterator corresponding to the start of the first group.
 * @param first_end    InputIterator corresponding to the end of the first group.
 * @param second_begin InputIterator corresponding to the start of the second group.
 * @param second_end   InputIterator corresponding to the end of the second group.
 *
 * @return Distance
 */
template <typename value_type, typename InputIterator>
value_type distance_correlation(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator second_end)
{
	GeneTrail::DenseMatrix A = pairwise_distance<value_type>(first_begin, first_end);
	A = normalize_distances(A);
	GeneTrail::DenseMatrix B = pairwise_distance<value_type>(second_begin, second_end);
	B = normalize_distances(B);
	double cov = distance_covariance<value_type>(A, B);
	double varA = distance_covariance<value_type>(A, A);
	double varB = distance_covariance<value_type>(B, B);
	return cov / std::sqrt(varA * varB);
}

/**
 * This methods implements the Dynamic Time Warping Distance Measure.
 *
 * @param first_begin  InputIterator corresponding to the start of the first group.
 * @param first_end    InputIterator corresponding to the end of the first group.
 * @param second_begin InputIterator corresponding to the start of the second group.
 * @param second_end   InputIterator corresponding to the end of the second group.
 *
 * @return Distance
 */
template <typename value_type, typename InputIterator>
value_type dynamic_time_warping(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator second_end)
{
	std::vector<double> x(first_begin, first_end);
	std::vector<double> y(second_begin, second_end);
	GeneTrail::DenseMatrix DTW(x.size(), y.size());
	DTW(0,0) = 0.0;
	for(size_t i=0; i<x.size(); ++i) {
                for(size_t j=0; j<y.size(); ++j) {
                        double cost = abs(x[i]-y[j]);
                        if(i == 0 && j == 0) continue;
                        if(i == 0) {
                                DTW(i,j) = cost + DTW(i, j-1);
                        } else if(j == 0) {
                                DTW(i,j) = cost + DTW(i-1, j);
                        } else {
                                DTW(i,j) = cost + std::min(DTW(i-1, j), std::min(DTW(i, j-1), DTW(i-1, j-1)));
                        }
                }
        }
        return DTW(x.size()-1, y.size()-1);
}

/**
 * This methods implements the Dynamic Time Warping Distance Measure (for gradients).
 *
 * @param first_begin  InputIterator corresponding to the start of the first group.
 * @param first_end    InputIterator corresponding to the end of the first group.
 * @param second_begin InputIterator corresponding to the start of the second group.
 * @param second_end   InputIterator corresponding to the end of the second group.
 *
 * @return Distance
 */
template <typename value_type, typename InputIterator>
value_type dynamic_time_warping_for_gradients(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator second_end)
{
	std::vector<double> x(first_begin, first_end);
	std::vector<double> y(second_begin, second_end);
	GeneTrail::DenseMatrix DTW(x.size()-1, y.size()-1);
	DTW(0,0) = 0.0;
	for(size_t i=0; i<x.size()-1; ++i) {
                for(size_t j=0; j<y.size()-1; ++j) {
                        double cost = abs((x[i+1]-x[i])-(y[j+1]-y[j]));
                        if(i == 0 && j == 0) continue;
                        if(i == 0) {
                                DTW(i,j) = cost + DTW(i, j-1);
                        } else if(j == 0) {
                                DTW(i,j) = cost + DTW(i-1, j);
                        } else {
                                DTW(i,j) = cost + std::min(DTW(i-1, j), std::min(DTW(i, j-1), DTW(i-1, j-1)));
                        }
                }
        }
        return DTW(x.size()-2, y.size()-2);
}
}
}

#endif // TRANSITIVITY_CLUSTERING_ILP_STATISTIC_H

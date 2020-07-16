/** 
 * Copyright (C) 2020 Tim Kehl <tkehl@bioinf.uni-sb.de>
 *                    Kerstin Lenhof <klenhof@bioinf.uni-sb.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public
 * License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
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
 * This methods implements the euclidean distance for gradients.
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
euclidean_distance_for_gradients(InputIterator first_begin, InputIterator first_end,
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
        return euclidean_distance<value_type>(fc_1.begin(), fc_1.end(), fc_2.begin(), fc_2.end());
}

template <typename value_type>
value_type theta(value_type x1, value_type x2, value_type y1, value_type y2)
{
	if(x1 == x2 && y1 == y2 && x1 == y1 && x2 == y2) {
		return 0.0;
	}
	value_type cos_alpha = (x1*y1 + x2*y2)/(std::sqrt(x1*x1 + x2*x2) * std::sqrt(y1*y1 + y2*y2));
	value_type theta = std::acos(cos_alpha) / 3.141593 * 180.0;
	return theta;
}


/**
 * This methods calculates the angles between ...
 *
 * @param first_begin  InputIterator corresponding to the start of the first group.
 * @param first_end    InputIterator corresponding to the end of the first group.
 * @param second_begin InputIterator corresponding to the start of the second group.
 * @param second_end   InputIterator corresponding to the end of the second group.
 *
 * @return Distance
 */
template <typename value_type, typename InputIterator>
std::vector<value_type> calculate_angles(InputIterator begin, InputIterator end)
{
	std::vector<value_type> angles;
	angles.reserve(std::distance(begin, end));
	for(; begin != (end-1); ++begin) {
		angles.emplace_back(theta(*(begin), *(begin), *(begin), *(begin+1)));
	}
	return std::move(angles);
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
angle_distance(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator)
{
	value_type sum = 0.0;
	for(auto i=0; i<std::distance(first_begin, first_end)-1; ++i) {
		double th = theta(1.0, *(first_begin+1+i) - *(first_begin+i), 1.0, *(second_begin+1+i) - *(second_begin+i));
		sum += th;
	}
	return sum / (double)std::distance(first_begin, first_end);
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
normalized_angle_distance(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator)
{
	value_type sum = 0.0;
	for(auto i=0; i<std::distance(first_begin, first_end)-1; ++i) {
		double fc1 = *(first_begin+1+i) - *(first_begin+i);
		double fc2 = *(second_begin+1+i) - *(second_begin+i);
		double th = theta(1.0, fc1, 1.0, fc2);
		double th1 = theta(1.0, 0.0, 1.0, fc1);
		double th2 = theta(1.0, 0.0, 1.0, fc2);
		sum += th / (1.0 + th1 + th2);
	}
	return sum / (double)std::distance(first_begin, first_end);
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
shifted_euclidean_distance_for_angles(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator second_end)
{
	std::vector<value_type> ang1 = calculate_angles<value_type>(first_begin, first_end);
	std::vector<value_type> ang2 = calculate_angles<value_type>(second_begin, second_end);
	return shifted_euclidean_distance_for_points<value_type>(ang1.begin(), ang1.end(), ang2.begin(), ang2.end());
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

/**
 * This methods calculates the pairwise euclidean distance for all points in a given range.
 *
 * @param begin  InputIterator corresponding to the start of the first group.
 * @param end    InputIterator corresponding to the end of the first group.
 *
 * @return Distance matrix
 */
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

/**
 * This methods normalizes a distance matrix.
 *
 * @param begin  InputIterator corresponding to the start of the first group.
 * @param end    InputIterator corresponding to the end of the first group.
 *
 * @return Normalized distance matrix
 */
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
 * @param A Normalized distance matrix
 * @param B Normalized distance matrix
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
                        double cost = std::abs(x[i]-y[j]);
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
                        double cost = std::abs((x[i+1]-x[i])-(y[j+1]-y[j]));
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

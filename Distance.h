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
                     InputIterator second_begin, InputIterator second_end)
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
                     InputIterator second_begin, InputIterator second_end)
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

}
}

#endif // TRANSITIVITY_CLUSTERING_ILP_STATISTIC_H

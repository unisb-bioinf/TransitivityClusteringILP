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

#ifndef GT2_CORE_STATISTIC_H
#define GT2_CORE_STATISTIC_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>
#include <iostream>
#include <map>

#include <boost/math/special_functions/atanh.hpp>

#include "macros.h"

namespace GeneTrail
{

/**
 * A collection of mathematical operations.
 */
namespace statistic
{
/**
 * This method calculates the absolute value for each entry of a given
 * range. The results are updated in place.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 */
template <typename value_type, typename InputIterator>
void abs(InputIterator begin, InputIterator end)
{
	std::transform(begin, end, begin,
	               static_cast<value_type (*)(value_type)>(std::abs));
}

/**
 * This method calculates the square root for each entry of a given
 * range. The results are updated in place.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 */
template <typename value_type, typename InputIterator>
void sqrt(InputIterator begin, InputIterator end)
{
	std::transform(begin, end, begin,
	               static_cast<value_type (*)(value_type)>(std::sqrt));
}

/**
 * This method calculates the natural logarithm for each entry of a
 * given range. The results are updated in place.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 */
template <typename value_type, typename InputIterator>
void log(InputIterator begin, InputIterator end)
{
	std::transform(begin, end, begin,
	               static_cast<value_type (*)(value_type)>(std::log));
}

/**
 * This method calculates the logarithm (base 10) for each entry of a
 * given range. The results are updated in place.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 */
template <typename value_type, typename InputIterator>
void log10(InputIterator begin, InputIterator end)
{
	std::transform(begin, end, begin,
	               static_cast<value_type (*)(value_type)>(std::log10));
}

/**
 * This method calculates the logarithm (base 2) for each entry of a
 * given range. The results are updated in place.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 */
template <typename value_type, typename InputIterator>
void log2(InputIterator begin, InputIterator end)
{
	std::transform(begin, end, begin,
	               static_cast<value_type (*)(value_type)>(std::log2));
}

/**
 * This method calculates the n-th power for each entry of a given
 * range. The results are updated in place.
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 */
template <typename value_type, typename InputIterator>
void pow(InputIterator begin, InputIterator end, int n)
{
	std::transform(begin, end, begin,
	               [n](value_type x) { return std::pow(x, n); });
}

/**
 * This method calculates the maximal value of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 *
 * @return Maximum of the given range
 */
template <typename value_type, typename InputIterator>
value_type max(InputIterator begin, InputIterator end)
{
	return *std::max_element(begin, end);
}

/**
 * This method calculates the minimal value of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 *
 * @return Minimum of the given range
 */
template <typename value_type, typename InputIterator>
value_type min(InputIterator begin, InputIterator end)
{
	return *std::min_element(begin, end);
}

/**
 * This method calculates the sum of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 *
 * @return Sum of the given range
 */
template <typename value_type, typename InputIterator>
value_type sum(InputIterator begin, InputIterator end)
{
	return std::accumulate(begin, end, 0.0);
}

/**
 * This method calculates the mean of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 *
 * @return Mean of the given range
 */
template <typename value_type, typename InputIterator>
double mean(InputIterator begin, InputIterator end)
{
	auto dist = std::distance(begin, end);
	if(dist == 0) {
		return value_type();
	}

	return (std::accumulate(begin, end, 0.0) / ((value_type)dist));
}

/**
 * This method calculates the absolute mean of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 *
 * @return Mean of the given range
 */
template <typename value_type, typename InputIterator>
double abs_mean(InputIterator begin, InputIterator end)
{
	auto dist = std::distance(begin, end);
	if(dist == 0) {
		return value_type();
	}

	std::vector<value_type> tmp(begin, end);
	abs<value_type>(tmp.begin(), tmp.end());

	return (std::accumulate(tmp.begin(), tmp.end(), 0.0) / ((value_type)dist));
}

/**
 * This method calculates the max-mean statistic of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 *
 * @return Max-mean statistic of the given range
 */
template <typename value_type, typename InputIterator>
value_type max_mean(InputIterator begin, InputIterator end)
{
	value_type n = std::distance(begin, end);

	if(n == 0) {
		return value_type();
	}

	value_type positive_sum = 0;
	value_type negative_sum = 0;
	for(auto iter = begin; iter != end; ++iter) {
		if(*iter < 0) {
			negative_sum -= *iter;
		} else {
			positive_sum += *iter;
		}
	}

	// In the original paper, they return absolute values.
	// But as we want to distinguish between positive and negative,
	// we return signed values.
	if(negative_sum > positive_sum) {
		return -negative_sum / n;
	} else {
		return positive_sum / n;
	}
}

/**
 * This method calculates the median of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 *
 * @return Median of the given range
 */
template <typename value_type, typename InputIterator>
value_type median(InputIterator begin, InputIterator end)
{
	std::vector<value_type> tmp(begin, end);

	const auto dist = std::distance(tmp.begin(), tmp.end());

	if(dist == 0) {
		return value_type();
	}

	const auto dist2 = dist / 2;
	auto median_position = tmp.begin() + dist2;
	std::nth_element(tmp.begin(), median_position, tmp.end());

	if(2 * dist2 == dist) {
		auto second = std::max_element(tmp.begin(), median_position);
		return (*median_position + *second) * value_type(0.5);
	} else {
		return *median_position;
	}
}

/**
 * This method returns the middle value of a given range.
 *
 * @param begin InputIterator
 * @param end InputIterator
 * @return Middel value of the given range
 */
template <typename value_type, typename InputIterator>
value_type middle(InputIterator begin, InputIterator end)
{
	std::vector<value_type> tmp(begin, end);

	const auto dist = std::distance(tmp.begin(), tmp.end());

	if(dist == 0) {
		return value_type();
	}

	const auto dist2 = dist / 2;
	auto median_position = tmp.begin() + dist2;
	std::nth_element(tmp.begin(), median_position, tmp.end());
	return *median_position;
}

/**
 * This method calculates the pooled variance of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples.
 * @param end   InputIterator corresponding to the end of the samples.
 *
 * @return Variance of the given range
 */
template <typename value_type, typename InputIterator>
value_type var(InputIterator begin, InputIterator end, value_type mean)
{
	int n = std::distance(begin, end);
	if(n <= 1) {
		return ((value_type)0);
	} else {
		value_type v = value_type();
		value_type ep = value_type();
		for(auto it = begin; it != end; ++it) {
			value_type diff = *it - mean;
			ep += diff;
			v += diff * diff;
		}
		return (v - ep * ep / n) / (n - 1);
	}
}

/**
 * This method calculates the pooled variance of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples
 * @param end   InputIterator corresponding to the end of the samples
 * @return Variance of the given range
 */
template <typename value_type, typename InputIterator>
value_type var(InputIterator begin, InputIterator end)
{
	int n = std::distance(begin, end);
	if(n <= 1) {
		return ((value_type)0);
	} else {
		const value_type s = std::accumulate(begin, end, value_type());
		const value_type mean = s / n;
		return var(begin, end, mean);
	}
}

/**
 * This implements the online algorithm by Knuth for calculating the variance.
 * The advantage of this method is, that every input value needs to be read only
 * once. It additionally computes the mean and the number of input samples.
 *
 * For more information see:
 * Donald Knuth, The Art of Computer Programming, Vol 2, page 232, 3rd edition.
 *
 * @param begin InputIterator pointing to the start of a range for which the
 *              variance should be computed.
 * @param end   InputIterator pointing at the end of the range.
 */
template <typename value_type, typename InputIterator>
std::tuple<value_type, value_type, size_t> mean_var_size(InputIterator begin,
                                                         InputIterator end)
{
	value_type mean = 0.0;
	value_type var = 0.0;
	size_t size = 0u;

	for(; begin != end; ++begin) {
		auto x = *begin;
		++size;
		auto delta = x - mean;
		mean += delta / size;
		var += delta * (x - mean);
	}

	if(size <= 1) {
		return std::make_tuple(mean, value_type(), size);
	}

	return std::make_tuple(mean, var / (size - 1), size);
}

/**
 * This method calculates the skewness of a given range.
 *
 * @param begin InputIterator
 * @param end InputIterator
 * @return Skewness of the given range
 */
template <typename value_type, typename InputIterator>
value_type skewness(InputIterator begin, InputIterator end)
{
	value_type mean = 0.0, var = 0.0;
	size_t size = 0u;
	std::tie(mean, var, size) = mean_var_size<value_type>(begin, end);
	value_type sd = std::sqrt(var);
	value_type skewness = 0;
	if(sd == 0) {
		return value_type();
	}

	for(; begin != end; ++begin) {
		skewness += std::pow((*begin - mean) / sd, 3);
	}

	return (skewness * size) / ((size - 1) * (size - 2));
}

/**
 * This method calculates the average absolute deviation of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples
 * @param end   InputIterator corresponding to the end of the samples
 * @param ave Average of the given range
 * @return Standard deviation of the given range
 */
template <typename value_type, typename InputIterator>
std::vector<value_type> absolute_deviations(InputIterator begin,
                                            InputIterator end, value_type ave)
{
	std::vector<value_type> ads;
	ads.reserve(std::distance(begin, end));

	for(auto it = begin; it != end; ++it) {
		ads.emplace_back(std::abs(*it - ave));
	}

	return ads;
}

/**
 * This method calculates the average absolute deviation of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples
 * @param end   InputIterator corresponding to the end of the samples
 * @return Standard deviation of the given range
 */
template <typename value_type, typename InputIterator>
value_type mean_absolute_deviation(InputIterator begin, InputIterator end)
{
	auto ads = absolute_deviations(begin, end, mean<value_type>(begin, end));
	return mean<value_type>(ads.begin(), ads.end());
}

/**
 * This method calculates the average absolute deviation of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples
 * @param end   InputIterator corresponding to the end of the samples
 * @return Standard deviation of the given range
 */
template <typename value_type, typename InputIterator>
value_type median_absolute_deviation(InputIterator begin, InputIterator end)
{
	auto ads = absolute_deviations(begin, end, median<value_type>(begin, end));
	return median<value_type>(ads.begin(), ads.end());
}

/**
 * This method calculates the standard deviation of a given range.
 *
 * @param begin InputIterator corresponding to the start of the samples
 * @param end   InputIterator corresponding to the end of the samples
 * @return Standard deviation of the given range
 */
template <typename value_type, typename InputIterator>
value_type sd(InputIterator begin, InputIterator end)
{
	return std::sqrt(var<value_type, InputIterator>(begin, end));
}

/**
 * This method calculates the sample covariance between two ranges of
 * values.
 *
 * @param first_begin  InputIterator corresponding to the start of the first
 *group.
 * @param first_end    InputIterator corresponding to the end of the first
 *group.
 * @param second_begin InputIterator corresponding to the start of the second
 *group.
 * @param second_end   InputIterator corresponding to the end of the second
 *group.
 *
 * @return Covariance of the given range
 */
template <typename value_type, typename InputIterator>
value_type cov(InputIterator first_begin, InputIterator first_end,
               InputIterator second_begin, InputIterator second_end,
               value_type mean1, value_type mean2)
{
	assert(std::distance(first_begin, first_end) ==
	       std::distance(second_begin, second_end));
	_unused(second_end);
	value_type cov = 0.0;

	const size_t n = std::distance(first_begin, first_end);
	for(; first_begin != first_end; ++first_begin, ++second_begin) {
		cov += (*first_begin - mean1) * (*second_begin - mean2);
	}

	return cov / (n - 1);
}

/**
 * This method calculates the sample covariance between two ranges of
 *values.
 *
 * @param first_begin InputIterator
 * @param first_end InputIterator
 * @param second_begin InputIterator
 * @param second_begin InputIterator
 * @return Covariance of the given range
 */
template <typename value_type, typename InputIterator>
value_type cov(InputIterator first_begin, InputIterator first_end,
               InputIterator second_begin, InputIterator second_end)
{
	assert(std::distance(first_begin, first_end) ==
	       std::distance(second_begin, second_end));
	value_type mean1 = mean<value_type, InputIterator>(first_begin, first_end);
	value_type mean2 =
	    mean<value_type, InputIterator>(second_begin, second_end);

	return cov(first_begin, first_end, second_begin, second_end, mean1, mean2);
}

/**
 * This methods implements Pearson's correlation coefficient.
 *
 * @param first_begin  InputIterator corresponding to the start of the first
 *group.
 * @param first_end    InputIterator corresponding to the end of the first
 *group.
 * @param second_begin InputIterator corresponding to the start of the second
 *group.
 * @param second_end   InputIterator corresponding to the end of the second
 *group.
 *
 * @return Correlation coefficient
 */
template <typename value_type, typename InputIterator>
value_type
pearson_correlation(InputIterator first_begin, InputIterator first_end,
                    InputIterator second_begin, InputIterator second_end)
{
	assert(std::distance(first_begin, first_end) ==
	       std::distance(second_begin, second_end));

	auto mvs1 = mean_var_size<value_type>(first_begin, first_end);
	auto mvs2 = mean_var_size<value_type>(second_begin, second_end);

	value_type covar = cov(first_begin, first_end, second_begin, second_end,
	                       std::get<0>(mvs1), std::get<0>(mvs2));

	return covar /
	       (std::sqrt(std::get<1>(mvs1)) * std::sqrt(std::get<1>(mvs2)));
}

/**
 */
template <typename value_type, typename InputIterator>
std::vector<int> ranks(InputIterator begin, InputIterator end)
{
	std::vector<value_type> tmp(begin, end);
	std::sort(tmp.begin(), tmp.end());
	std::vector<int> ranks;
	for(auto it = begin; it != end; ++it) {
		for(size_t i = 0; i < tmp.size(); ++i) {
			if(*it == *(tmp.begin() + i)) {
				ranks.push_back(i);
			}
		}
	}
	return ranks;
}

/**
 * This methods implements Spearman's correlation coefficient.
 *
 * @param first_begin  InputIterator corresponding to the start of the first
 *group.
 * @param first_end    InputIterator corresponding to the end of the first
 *group.
 * @param second_begin InputIterator corresponding to the start of the second
 *group.
 * @param second_end   InputIterator corresponding to the end of the second
 *group.
 *
 * @return Correlation coefficient
 */
template <typename value_type, typename InputIterator>
value_type
spearman_correlation(InputIterator first_begin, InputIterator first_end,
                     InputIterator second_begin, InputIterator second_end)
{
	std::vector<int> first_ranks =
	    ranks<value_type, InputIterator>(first_begin, first_end);
	std::vector<int> second_ranks =
	    ranks<value_type, InputIterator>(second_begin, second_end);
	return pearson_correlation<value_type, std::vector<int>::iterator>(
	    first_ranks.begin(), first_ranks.end(), second_ranks.begin(),
	    second_ranks.end());
}

/**
 * This method implements the Kendall Tau correlation coefficient.
 *
 * @param first_begin InputIterator
 * @param first_end InputIterator
 * @param second_begin InputIterator
 * @param second_begin InputIterator
 * @return correlation coefficient
 */
template <typename value_type, typename InputIterator>
value_type
kendall_tau_correlation(InputIterator first_begin, InputIterator first_end,
                        InputIterator second_begin, InputIterator second_end)
{
	int dist_second = std::distance(second_begin, second_end);

	assert(std::distance(first_begin, first_end) == dist_second);

	std::vector<int> first_ranks =
	    ranks<value_type, InputIterator>(first_begin, first_end);

	std::map<int, int> rank_map;

	int counter = 0;
	for(int rank : first_ranks) {
		rank_map[rank] = counter++;
	}

	int greater_count = 0, less_count = 0, x_equal_count = 0, y_equal_count = 0;

	for(int i = 0; i < dist_second; i++) {
		for(int j = i + 1; j < dist_second; j++) {
			auto x_i = *(first_begin + rank_map[i]);
			auto x_j = *(first_begin + rank_map[j]);
			auto y_i = *(second_begin + rank_map[i]);
			auto y_j = *(second_begin + rank_map[j]);

			if(x_i == x_j) {
				if(y_i != y_j) {
					++x_equal_count;
				}
				continue;
			}

			if(y_i == y_j) {
				++y_equal_count;
				continue;
			}

			if(y_i < y_j) {
				++less_count;
				continue;
			}

			++greater_count;
		}
	}

	value_type denominator =
	    std::sqrt((greater_count + less_count + x_equal_count) *
	              (greater_count + less_count + y_equal_count));

	assert(denominator != 0);

	return (less_count - greater_count) / denominator;
}

/**
 * This function converts the given correlation coefficient
 * into a t-score with n-2 degrees of freedom.
 *
 * @param r Correlation coefficient
 * @param n Sample size
 *
 * @return A t-score that follows a t-distribution with n-2 degrees of freedom.
 */
  template <typename value_type>
  value_type correlation_to_t(value_type r, value_type n) {
	  return r / ((1.0 - r*r) / (n-2));
  }

/**
 * This function implements the Fisher's z-transformation 
 * for correlation coefficients.
 *
 * @param r Correlation coefficient
 *
 * @return Z-score that follows a normal distribution with sd 1 / sqrt(n-3).
 */
  template <typename value_type>
  value_type fisher_z_transform(value_type r) {
	  return boost::math::atanh(r);
  }

/**
 * This function calculates the standard error for 
 * a Fisher's z-transformed correlation coefficient
 * with n samples.
 *
 * @param n Sample size
 *
 * @return Standard error.
 */
  template <typename value_type>
  value_type fisher_z_transform_sd(value_type n) {
	  return 1.0 / sqrt(n - 3.0);
  }

/**
 * This function converts a correlation coefficient
 * into a z-score using Fisher's z-transformation.
 *
 * @param n Sample size
 *
 * @return Standard error.
 */
  template <typename value_type>
  value_type correlation_to_z(value_type r, value_type n, value_type mean = 0.0) {
	  return (fisher_z_transform(r) - mean) / fisher_z_transform_sd(n);
  }

/**
 * This function calculates log(mean(a)) -  log(mean(b)).
 * Log-Mean-Fold-Quotient
 *
 * @param begin1 InputIterator corresponding to the start of the first group.
 * @param end1   InputIterator corresponding to the end of the first group.
 * @param begin2 InputIterator corresponding to the start of the second group.
 * @param end2   InputIterator corresponding to the end of the second group.
 *
 * @return log(mean(a)) -  log(mean(b))
 */
template <typename value_type, typename InputIterator1, typename InputIterator2>
value_type log_mean_fold_quotient(InputIterator1 begin1, InputIterator1 end1,
                                  InputIterator2 begin2, InputIterator2 end2)
{
	return std::log(mean<value_type, InputIterator1>(begin1, end1)) -
	       std::log(mean<value_type, InputIterator2>(begin2, end2));
}

/**
 * Log-Median-Fold-Quotient
 * This function calculates log(median(a)) -  log(median(b)).
 *
 * @param begin1 InputIterator corresponding to the start of the first group.
 * @param end1   InputIterator corresponding to the end of the first group.
 * @param begin2 InputIterator corresponding to the start of the second group.
 * @param end2   InputIterator corresponding to the end of the second group.
 *
 * @return log(median(a)) -  log(median(b))
 */
template <typename value_type, typename InputIterator1, typename InputIterator2>
value_type log_median_fold_quotient(InputIterator1 begin1, InputIterator1 end1,
                                    InputIterator2 begin2, InputIterator2 end2)
{
	return std::log(median<value_type, InputIterator1>(begin1, end1)) -
	       std::log(median<value_type, InputIterator2>(begin2, end2));
}

/**
 * This function calculates mean(a) /  mean(b).
 * Mean-Fold-Quotient
 *
 * @param begin1 InputIterator corresponding to the start of the first group.
 * @param end1   InputIterator corresponding to the end of the first group.
 * @param begin2 InputIterator corresponding to the start of the second group.
 * @param end2   InputIterator corresponding to the end of the second group.
 *
 * @return mean(a) /  mean(b)
 */
template <typename value_type, typename InputIterator1, typename InputIterator2>
value_type mean_fold_quotient(InputIterator1 begin1, InputIterator1 end1,
                              InputIterator2 begin2, InputIterator2 end2)
{
	return mean<value_type, InputIterator1>(begin1, end1) /
	       mean<value_type, InputIterator2>(begin2, end2);
}

/**
* This function calculates mean(a) -  mean(b).
* Mean-Fold-Difference
*
 * @param begin1 InputIterator corresponding to the start of the first group.
 * @param end1   InputIterator corresponding to the end of the first group.
 * @param begin2 InputIterator corresponding to the start of the second group.
 * @param end2   InputIterator corresponding to the end of the second group.
 *
* @return mean(a) -  mean(b)
*/
template <typename value_type, typename InputIterator1, typename InputIterator2>
value_type mean_fold_difference(InputIterator1 begin1, InputIterator1 end1,
                                InputIterator2 begin2, InputIterator2 end2)
{
	return mean<value_type, InputIterator1>(begin1, end1) -
	       mean<value_type, InputIterator2>(begin2, end2);
}

/**
 * Median-Fold-Quotient
 * This function calculates median(a) /  median(b).
 *
 * @param begin1 InputIterator corresponding to the start of the first group.
 * @param end1   InputIterator corresponding to the end of the first group.
 * @param begin2 InputIterator corresponding to the start of the second group.
 * @param end2   InputIterator corresponding to the end of the second group.
 *
 * @return median(a) /  median(b)
 */
template <typename value_type, typename InputIterator1, typename InputIterator2>
value_type median_fold_quotient(InputIterator1 begin1, InputIterator1 end1,
                                InputIterator2 begin2, InputIterator2 end2)
{
	return median<value_type, InputIterator1>(begin1, end1) /
	       median<value_type, InputIterator2>(begin2, end2);
}

/**
 * Mean-Fold-Difference
 * This function calculates median(a) -  median(b).
 *
 * @param begin1 InputIterator corresponding to the start of the first group.
 * @param end1   InputIterator corresponding to the end of the first group.
 * @param begin2 InputIterator corresponding to the start of the second group.
 * @param end2   InputIterator corresponding to the end of the second group.
 *
 * @return median(a) -  median(b)
 */
template <typename value_type, typename InputIterator1, typename InputIterator2>
value_type median_fold_difference(InputIterator1 begin1, InputIterator1 end1,
                                  InputIterator2 begin2, InputIterator2 end2)
{
	return median<value_type, InputIterator1>(begin1, end1) -
	       median<value_type, InputIterator2>(begin2, end2);
}

/**
 * This function calculates the z-score.
 *
 * @param x Input for which the z-score is requested
 * @param begin InputIterator corresponding to the start of the samples
 * @param end   InputIterator corresponding to the end of the samples
 *
 * @return Z-score
 */
template <typename value_type, typename InputIterator1>
value_type z_score(value_type x, InputIterator1 begin, InputIterator1 end)
{
	return (x - mean<value_type, InputIterator1>(begin, end)) /
	       sd<value_type, InputIterator1>(begin, end);
}
}
}

#endif // GT2_CORE_STATISTIC_H

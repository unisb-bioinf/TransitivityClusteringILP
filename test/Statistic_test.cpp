#include <gtest/gtest.h>

#include "../Statistic.h"

#include <iostream>
#include <initializer_list>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace GeneTrail;

double TOLERANCE = 0.00001;


constexpr std::initializer_list<double> ai = {35.5,-31.7,31.2,36.6,-22.8,28.0,-24.6,26.1,-34.5,27.7};
constexpr std::initializer_list<double> bi = {45.3,-36.0,38.6,-44.7,31.4,-33.5,28.8,-35.8,42.9,-35.0};
constexpr std::initializer_list<double> ci = {24.5,-23.5,12.5,-2.0,10.0,5.5,4.5,-3.5,-21.5,9.0,-11.5,1.5,-1.5,4.0,12.0};

std::vector<double> a(ai);
std::vector<double> b(bi);
std::vector<double> c(ci);

TEST(Statistic, Max)
{
	double max1 = statistic::max<double,std::vector<double>::iterator>(a.begin(),a.end());
	EXPECT_NEAR(max1, 36.6, TOLERANCE);
	double max2 = statistic::max<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(max2, 45.3, TOLERANCE);
}

TEST(Statistic, Min)
{
	double min1 = statistic::min<double, std::vector<double>::iterator>(a.begin(),a.end());
	EXPECT_NEAR(min1, -34.5, TOLERANCE);
	double min2 = statistic::min<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(min2, -44.7, TOLERANCE);
}

TEST(Statistic, Abs)
{
	auto tmp = a;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	for(size_t i = 0; i < a.size(); ++i) {
		EXPECT_EQ(tmp[i],std::abs(a[i]));
	}
	tmp = b;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	for(size_t i = 0; i < b.size(); ++i) {
		EXPECT_EQ(tmp[i],std::abs(b[i]));
	}
}

TEST(Statistic, Mean)
{
	double mean1 = statistic::mean<double, std::vector<double>::iterator>(a.begin(),a.end());
	EXPECT_NEAR(mean1, 7.15, TOLERANCE);
	double mean2 = statistic::mean<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(mean2, 0.2, TOLERANCE);
}

TEST(Statistic, Median_Length_0)
{
	std::vector<double> values;
	double median = statistic::median<double>(values.begin(), values.end());
	EXPECT_DOUBLE_EQ(0.0, median);
}


TEST(Statistic, Median_Length_1)
{
	std::vector<double> values { 1.0 };
	double median = statistic::median<double>(values.begin(), values.end());
	EXPECT_DOUBLE_EQ(1.0, median);
}

TEST(Statistic, Median_Length_2)
{
	std::vector<double> values { -1.0, 1.0 };
	double median = statistic::median<double>(values.begin(), values.end());
	EXPECT_DOUBLE_EQ(0.0, median);
	std::vector<double> values2 { 1.0, 3.0 };
	double median2 = statistic::median<double>(values2.begin(), values2.end());
	EXPECT_DOUBLE_EQ(2.0, median2);
}

TEST(Statistic, Median)
{
	double median1 = statistic::median<double>(a.begin(),a.end());
	EXPECT_NEAR(median1, 26.9, TOLERANCE);
	double median2 = statistic::median<double>(b.begin(),b.end());
	EXPECT_NEAR(median2, -2.35, TOLERANCE);
	double median3 = statistic::median<double>(c.begin(),c.end());
	EXPECT_DOUBLE_EQ(4.0, median3); // Odd length vectors should be exact
}

TEST(Statistic, Middle)
{
	double median1 = statistic::middle<double>(a.begin(),a.end());
	EXPECT_EQ(median1, 27.7);
	double median2 = statistic::middle<double>(b.begin(),b.end());
	EXPECT_EQ(median2, 28.8);
	double median3 = statistic::middle<double>(c.begin(),c.end());
	EXPECT_EQ(4, median3);
}

TEST(Statistic, Var)
{
	double var1 = statistic::var<double, std::vector<double>::iterator>(a.begin(),a.end());
	EXPECT_NEAR(var1, 957.185, TOLERANCE);
	double var2 = statistic::var<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(var2, 1568.93778, TOLERANCE);
}

TEST(Statistic, Sd)
{
	double sd1 = statistic::sd<double, std::vector<double>::iterator>(a.begin(),a.end());
	EXPECT_NEAR(sd1, 30.93841, TOLERANCE);
	double sd2 = statistic::sd<double, std::vector<double>::iterator>(b.begin(),b.end());
	EXPECT_NEAR(sd2, 39.60982, TOLERANCE);
}

TEST(Statistic, skewness)
{
	double skew1 = statistic::skewness<double, std::vector<double>::iterator>(a.begin(), a.end());
	EXPECT_NEAR(skew1, -0.4829496, TOLERANCE);
	double skew2 = statistic::skewness<double, std::vector<double>::iterator>(b.begin(), b.end());
	EXPECT_NEAR(skew2, 0.0301191, TOLERANCE);
	double skew3 = statistic::skewness<double, std::vector<double>::iterator>(c.begin(), c.end());
	EXPECT_NEAR(skew3, -0.5220547, TOLERANCE);
}

TEST(Statistic, MeanAbsoluteDeviation)
{
	double mad1 = statistic::mean_absolute_deviation<double, std::vector<double>::iterator>(a.begin(), a.end());
	EXPECT_NEAR(mad1, 28.44, TOLERANCE);
	double mad2 = statistic::mean_absolute_deviation<double, std::vector<double>::iterator>(b.begin(), b.end());
	EXPECT_NEAR(mad2, 37.2, TOLERANCE);
	double mad3 = statistic::mean_absolute_deviation<double, std::vector<double>::iterator>(c.begin(), c.end());
	EXPECT_NEAR(mad3, 9.533333, TOLERANCE);
}

TEST(Statistic, MedianAbsoluteDeviation)
{
	double mad1 = statistic::median_absolute_deviation<double, std::vector<double>::iterator>(a.begin(), a.end());
	EXPECT_NEAR(mad1, 9.15, TOLERANCE);
	double mad2 = statistic::median_absolute_deviation<double, std::vector<double>::iterator>(b.begin(), b.end());
	EXPECT_NEAR(mad2, 33.7, TOLERANCE);
	double mad3 = statistic::median_absolute_deviation<double, std::vector<double>::iterator>(c.begin(), c.end());
	EXPECT_NEAR(mad3, 6, TOLERANCE);
}

TEST(Statistic, Log){
	auto tmp = a;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	statistic::log<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	for(size_t i = 0; i < a.size(); ++i) {
		EXPECT_EQ(tmp[i], std::log(std::abs(a[i])));
	}
	tmp = b;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	statistic::log<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	for(size_t i = 0; i < b.size(); ++i) {
		EXPECT_EQ(tmp[i], std::log(std::abs(b[i])));
	}
}

TEST(Statistic, Log2)
{
	auto tmp = a;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	statistic::log2<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	for(size_t i = 0; i < a.size(); ++i) {
		EXPECT_EQ(tmp[i], std::log2(std::abs(a[i])));
	}
	tmp = b;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	statistic::log2<double, std::vector<double>::iterator>(tmp.begin(),tmp.end());
	for(size_t i = 0; i < b.size(); ++i) {
		EXPECT_EQ(tmp[i], std::log2(std::abs(b[i])));
	}
}

TEST(Statistic, Log10)
{
	auto tmp = a;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),
	                                                      tmp.end());
	statistic::log10<double, std::vector<double>::iterator>(tmp.begin(),
	                                                       tmp.end());
	for(size_t i = 0; i < a.size(); ++i) {
		EXPECT_EQ(tmp[i], std::log10(std::abs(a[i])));
	}
	tmp = b;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),
	                                                      tmp.end());
	statistic::log10<double, std::vector<double>::iterator>(tmp.begin(),
	                                                       tmp.end());
	for(size_t i = 0; i < b.size(); ++i) {
		EXPECT_EQ(tmp[i], std::log10(std::abs(b[i])));
	}
}

TEST(Statistic, Sqrt)
{
	auto tmp = a;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),
	                                                      tmp.end());
	statistic::sqrt<double, std::vector<double>::iterator>(tmp.begin(),
	                                                       tmp.end());
	for(size_t i = 0; i < a.size(); ++i) {
		EXPECT_EQ(tmp[i], std::sqrt(std::abs(a[i])));
	}
	tmp = b;
	statistic::abs<double, std::vector<double>::iterator>(tmp.begin(),
	                                                      tmp.end());
	statistic::sqrt<double, std::vector<double>::iterator>(tmp.begin(),
	                                                       tmp.end());
	for(size_t i = 0; i < b.size(); ++i) {
		EXPECT_EQ(tmp[i], std::sqrt(std::abs(b[i])));
	}
}

TEST(Statistic, Pow)
{
	auto tmp = a;
	statistic::pow<double, std::vector<double>::iterator>(tmp.begin(), tmp.end(),2);
	for(size_t i = 0; i < a.size(); ++i) {
		EXPECT_EQ(tmp[i], a[i] * a[i]);
	}
	tmp = b;
	statistic::pow<double, std::vector<double>::iterator>(tmp.begin(), tmp.end(),2);
	for(size_t i = 0; i < b.size(); ++i) {
		EXPECT_EQ(tmp[i], b[i] * b[i]);
	}
}

TEST(Statistic, Cov)
{
	auto tmp = a;
	auto tmp2 = b;
	auto covar = statistic::cov<double, std::vector<double>::iterator>(tmp.begin(), tmp.end(), tmp2.begin(), tmp2.end());
	EXPECT_NEAR(covar, -382.553333, TOLERANCE);
}

TEST(Statistic, Pearson)
{
	auto tmp = a;
	auto tmp2 = b;
	auto covar = statistic::pearson_correlation<double, std::vector<double>::iterator>(tmp.begin(), tmp.end(), tmp2.begin(), tmp2.end());
	EXPECT_NEAR(covar, -0.31217, TOLERANCE);
}

TEST(Statistic, Spearman)
{
	auto tmp = a;
	auto tmp2 = b;
	auto covar = statistic::spearman_correlation<double, std::vector<double>::iterator>(tmp.begin(), tmp.end(), tmp2.begin(), tmp2.end());
	EXPECT_NEAR(covar, -0.06666667, TOLERANCE);
}

TEST(Statistic, Kendall)
{
	std::initializer_list<double> a = {9.694027318, 9.597726177, 9.179537375, 9.027478569, 8.792238023, 7.40578791, 7.104261273, 8.594364081, 7.302828273, 7.288458962};
	std::initializer_list<double> b1 = {5.747209717, 5.819584082, 5.616658389, 5.577300854, 5.559250186, 5.516019286, 5.41052816, 5.561140502, 5.567504767, 5.545473661};
	std::initializer_list<double> b2 = {7.793992658, 7.518006447, 7.41210759, 7.260149192, 7.511888682, 7.749760069, 7.502800333, 7.672162745, 7.971909606, 7.621752909};
	
	std::vector<double> a_v(a);
	std::vector<double> b1_v(b1);
	std::vector<double> b2_v(b2);
	
	auto kendall1 = statistic::kendall_tau_correlation<double, std::vector<double>::iterator>(a_v.begin(), a_v.end(), b1_v.begin(), b1_v.end());
	EXPECT_NEAR(kendall1, 0.73333, TOLERANCE);
	
	auto kendall2 = statistic::kendall_tau_correlation<double, std::vector<double>::iterator>(a_v.begin(), a_v.end(), b2_v.begin(), b2_v.end());
	EXPECT_NEAR(kendall2, -0.06666, TOLERANCE);
}

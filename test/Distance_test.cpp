#include <gtest/gtest.h>

#include "../Distance.h"

#include <iostream>
#include <initializer_list>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace TransitivityClusteringILP::Distance;

double TOLERANCE = 0.00001;


constexpr std::initializer_list<double> ai = {35.5,-31.7,31.2,36.6,-22.8,28.0,-24.6,26.1,-34.5,27.7};
constexpr std::initializer_list<double> bi = {45.3,-36.0,38.6,-44.7,31.4,-33.5,28.8,-35.8,42.9,-35.0};
constexpr std::initializer_list<double> ci = {24.5,-23.5,12.5,-2.0,10.0,5.5,4.5,-3.5,-21.5,9.0};

std::vector<double> a(ai);
std::vector<double> b(bi);
std::vector<double> c(ci);

TEST(Statistic, SquaredDistance)
{
	double d1 = squared_distance<double>(a.begin(), a.end(), b.begin(), b.end());
	EXPECT_NEAR(d1, 3010.409, TOLERANCE);
	double d2 = squared_distance<double>(a.begin(), a.end(), c.begin(), c.end());
        EXPECT_NEAR(d2, 585.164, TOLERANCE);
	double d3 = squared_distance<double>(c.begin(), c.end(), b.begin(), b.end());
        EXPECT_NEAR(d3, 1278.949, TOLERANCE);
}

/**
eudist <- function(x, y){
  diff <- x - y
  dist <- sqrt(mean(diff * diff))
  return(dist)
}
*/
TEST(Statistic, EuclideanDistance)
{
        double d1 = euclidean_distance<double>(a.begin(), a.end(), b.begin(), b.end());
        EXPECT_NEAR(d1, 54.86719, TOLERANCE);
        double d2 = euclidean_distance<double>(a.begin(), a.end(), c.begin(), c.end());
        EXPECT_NEAR(d2, 24.19016, TOLERANCE);
        double d3 = euclidean_distance<double>(c.begin(), c.end(), b.begin(), b.end());
        EXPECT_NEAR(d3, 35.7624, TOLERANCE);
}


/*

seudist <- function(x, y){
    diff <- x - y
    s <- mean(diff)
    diff <- x - s - y
    dist <- sqrt(mean(diff * diff))
    return(dist)
}

*/

TEST(Statistic, ShiftedEuclideanDistanceForPoints)
{
        double d1 = shifted_euclidean_distance_for_points<double>(a.begin(), a.end(), b.begin(), b.end());
        EXPECT_NEAR(d1, 54.42524, TOLERANCE);
        double d2 = shifted_euclidean_distance_for_points<double>(a.begin(), a.end(), c.begin(), c.end());
        EXPECT_NEAR(d2, 23.53304, TOLERANCE);
        double d3 = shifted_euclidean_distance_for_points<double>(c.begin(), c.end(), b.begin(), b.end());
        EXPECT_NEAR(d3, 35.73691, TOLERANCE);
}

/*

fcs <- function(x) {
    g <- c()
    for(i in (c(1:length(x)-1))){
        g <- c(g, x[i+1] - x[i])
    }
    return(g)
}

*/

TEST(Statistic, ShiftedEuclideanDistanceForGradients)
{
        double d1 = shifted_euclidean_distance_for_gradients<double>(a.begin(), a.end(), b.begin(), b.end());
        EXPECT_NEAR(d1, 107.9582, TOLERANCE*10);
        double d2 = shifted_euclidean_distance_for_gradients<double>(a.begin(), a.end(), c.begin(), c.end());
        EXPECT_NEAR(d2, 45.41857, TOLERANCE);
        double d3 = shifted_euclidean_distance_for_gradients<double>(c.begin(), c.end(), b.begin(), b.end());
        EXPECT_NEAR(d3, 69.06889, TOLERANCE);
}

TEST(Statistic, DistanceCovariance)
{
	double x = distance_covariance<double>(a.begin(), a.end(), a.begin(), a.end());
	EXPECT_NEAR(x, 26.39902, TOLERANCE);
        double d1 = distance_covariance<double>(a.begin(), a.end(), b.begin(), b.end());
        EXPECT_NEAR(d1, 11.89224, TOLERANCE);
        double d2 = distance_covariance<double>(a.begin(), a.end(), c.begin(), c.end());
        EXPECT_NEAR(d2, 9.963814, TOLERANCE);
        double d3 = distance_covariance<double>(c.begin(), c.end(), b.begin(), b.end());
        EXPECT_NEAR(d3, 7.9316, TOLERANCE);
}

/*TEST(Statistic, DistanceCorrelation)
{
        double d1 = distance_correlation<double>(a.begin(), a.end(), b.begin(), b.end());
        EXPECT_NEAR(d1, 0.3923192, TOLERANCE);
        double d2 = distance_correlation<double>(a.begin(), a.end(), c.begin(), c.end());
        EXPECT_NEAR(d2, 0.6148018, TOLERANCE);
        double d3 = distance_correlation<double>(c.begin(), c.end(), b.begin(), b.end());
        EXPECT_NEAR(d3, 0.4262202, TOLERANCE);
}*/

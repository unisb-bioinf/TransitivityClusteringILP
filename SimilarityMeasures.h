#ifndef SIMILARITY_MEASURES_H
#define SIMILARITY_MEASURES_H

#include "Statistic.h"
#include "Distance.h"

namespace TransitivityClusteringILP {

class PearsonCorrelation
{
  public:
	template <typename Iterator>
	double compute_similarity(Iterator first_begin, Iterator first_end, Iterator second_begin, Iterator second_end)
	{
		return GeneTrail::statistic::pearson_correlation<double>(first_begin, first_end, second_begin, second_end);
	}
};

class SpearmanCorrelation
{ 
  public:
        template <typename Iterator>
        double compute_similarity(Iterator first_begin, Iterator first_end, Iterator second_begin, Iterator second_end)
        {       
                return GeneTrail::statistic::spearman_correlation<double>(first_begin, first_end, second_begin, second_end);
        }
};

class DistanceCorrelation
{ 
  public:
        template <typename Iterator>
        double compute_similarity(Iterator first_begin, Iterator first_end, Iterator second_begin, Iterator second_end)
        {       
                return Distance::distance_correlation<double>(first_begin, first_end, second_begin, second_end);
        }
};

class EuclideanDistance
{ 
  public:
        template <typename Iterator>
        double compute_similarity(Iterator first_begin, Iterator first_end, Iterator second_begin, Iterator second_end)
        {       
                double dist = Distance::euclidean_distance<double>(first_begin, first_end, second_begin, second_end);
		return 1 / (1 + dist);
        }
};

class ShiftedEuclideanDistanceForPoints
{ 
  public:
        template <typename Iterator>
        double compute_similarity(Iterator first_begin, Iterator first_end, Iterator second_begin, Iterator second_end)
        {       
                double dist = Distance::shifted_euclidean_distance_for_points<double>(first_begin, first_end, second_begin, second_end);
		return 1 / (1 + dist);
        }
};

class ShiftedEuclideanDistanceForGradients
{ 
  public:
        template <typename Iterator>
        double compute_similarity(Iterator first_begin, Iterator first_end, Iterator second_begin, Iterator second_end)
        {       
                double dist = Distance::shifted_euclidean_distance_for_gradients<double>(first_begin, first_end, second_begin, second_end);
		return 1 / (1 + dist);
        }
};

class ShiftedEuclideanDistanceForGradientsAndPoints
{ 
  public:
        template <typename Iterator>
        double compute_similarity(Iterator first_begin, Iterator first_end, Iterator second_begin, Iterator second_end)
        {       
                double dist = Distance::shifted_euclidean_distance_for_gradients_and_points<double>(first_begin, first_end, second_begin, second_end);
		return 1 / (1 + dist);
        }
};

class DynamicTimeWarping
{ 
  public:
        template <typename Iterator>
        double compute_similarity(Iterator first_begin, Iterator first_end, Iterator second_begin, Iterator second_end)
        {       
                double dist = Distance::dynamic_time_warping<double>(first_begin, first_end, second_begin, second_end);
		return 1 / (1 + dist);
        }
};

class DynamicTimeWarpingForGradients
{ 
  public:
        template <typename Iterator>
        double compute_similarity(Iterator first_begin, Iterator first_end, Iterator second_begin, Iterator second_end)
        {       
                double dist = Distance::dynamic_time_warping_for_gradients<double>(first_begin, first_end, second_begin, second_end);
		return 1 / (1 + dist);
        }
};

}

#endif // SIMILARITY_MEASURES_H

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <string>

#include <boost/program_options.hpp>

#include "DenseMatrixReader.h"
#include "DenseMatrix.h"
#include "matrixTools.h"
#include "MatrixIterator.h"
#include "SimilarityMeasures.h"
#include "ClusteringILP2.h"
#include "ClusterEvaluation.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string matrix_ = "", similarity_measure_ = "";
bool isDistanceMatrix;
double similarity_threshold_;

MatrixReaderOptions matrixOptions;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("matrix,m", bpo::value<std::string>(&matrix_)->required(), "Name of a text file containing expression values as a matrix.")
		("distance,d", bpo::value<bool>(&isDistanceMatrix)->default_value(false)->zero_tokens(), "The matrix already is a distance matrix.")
		("no-row-names,r", bpo::value<bool>(&matrixOptions.no_rownames)->default_value(false)->zero_tokens(), "Does the file contain row names.")
		("no-col-names,c", bpo::value<bool>(&matrixOptions.no_colnames)->default_value(false)->zero_tokens(), "Does the file contain column names.")
		("add-col-name,a", bpo::value<bool>(&matrixOptions.additional_colname)->default_value(false)->zero_tokens(), "File containing two lines specifying which rownames belong to which group.")
		("similarity-measure,s", bpo::value<std::string>(&similarity_measure_)->default_value("spearman-correlation"), "Method used to compute pairwise similarities.")
		("similarity-threshold,t", bpo::value<double>(&similarity_threshold_)->default_value(0.95), "Threshold.");
	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
		bpo::notify(vm);
	}
	catch(bpo::error& e)
	{
		std::cerr << "ERROR: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return true;
}

template<typename SimilarityMeasure>
DenseMatrix compute_similarity_matrix(const DenseMatrix& matrix, SimilarityMeasure meas) {
	DenseMatrix dist(matrix.rows(), matrix.rows());
	dist.setRowNames(matrix.rowNames());
	dist.setColNames(matrix.rowNames());
	for(size_t i=0; i<matrix.rows(); ++i){
		dist.set(i, i, 1.0);
		for(size_t j=i+1; j<matrix.rows(); ++j){
			RowMajorMatrixIterator<Matrix> iit(&matrix, i), jit(&matrix, j);
			double d = meas.compute_similarity(
				iit->begin(),
				iit->end(),
				jit->begin(),
				jit->end()
			);
			dist.set(i, j, d);
			dist.set(j, i, d);
		}
	}
	return std::move(dist);
}

DenseMatrix compute_similarity_matrix(const DenseMatrix& matrix, const std::string& similatity_measure) {
	if(similatity_measure == "pearson-correlation") {
		return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::PearsonCorrelation()));
	} else if (similatity_measure == "spearman-correlation") {
                return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::SpearmanCorrelation()));
        } else if (similatity_measure == "distance-correlation") {
		return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::DistanceCorrelation()));
	} else if (similatity_measure == "euclidean-distance") {
                return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::EuclideanDistance()));
        } else if (similatity_measure == "shifted-euclidean-distance-for-points") {
		return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::ShiftedEuclideanDistanceForPoints()));
	} else if (similatity_measure == "shifted-euclidean-distance-for-gradients") {
                return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::ShiftedEuclideanDistanceForGradients()));
        } else if (similatity_measure == "shifted-euclidean-distance-for-gradients-and-points"){
                return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::ShiftedEuclideanDistanceForGradientsAndPoints()));
        } else if (similatity_measure == "dynamic-time-warping"){
                return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::DynamicTimeWarping()));
        } else if (similatity_measure == "dynamic-time-warping-for-gradients"){
                return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::DynamicTimeWarpingForGradients()));
        } else {
		return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::PearsonCorrelation()));
	}
}

std::vector<std::tuple<size_t, size_t, double>> extract_sorted_edges(const DenseMatrix& matrix) {
	std::vector<std::tuple<size_t, size_t, double>> edges;
	edges.reserve(matrix.rows()*matrix.rows()/2 - matrix.rows());
	for(size_t i=0; i<matrix.rows(); ++i){
                for(size_t j=i+1; j<matrix.rows(); ++j){
			edges.emplace_back(std::make_tuple(i, j, matrix(i,j)));
		}
	}
	std::sort(edges.begin(), edges.end(), [](const auto& a, const auto& b){return std::get<2>(a) > std::get<2>(b);});
	return std::move(edges);
}

int main(int argc, char* argv[])
{
	if(!parseArguments(argc, argv))
	{
		return -1;
	}

	if(matrixOptions.additional_colname && matrixOptions.no_rownames) {
		std::cerr << "Conflicting arguments. Additional colnames can only be "
		             "specified if row names are present!" << std::endl;
		return -2;
	}

	std::cout << "INFO: Reading matrix ..." << std::endl;
	DenseMatrix matrix(0,0);
	DenseMatrix sim(0,0);
	try {
		matrix = readDenseMatrix(matrix_, matrixOptions);
	} catch(const IOError& e) {
		std::cerr << "ERROR: Could not open input data matrix for reading." << std::endl;
		return -4;
	}

	
	if(!isDistanceMatrix) {
		std::cout << "INFO: Calculating similatiry matrix ..." << std::endl;
		sim = compute_similarity_matrix(matrix, similarity_measure_);
	} else {
		sim = readDenseMatrix(matrix_, matrixOptions);
        } 

	std::cout << "INFO: Extractin sorted edges ..." << std::endl;
	std::vector<std::tuple<size_t, size_t, double>> edges = extract_sorted_edges(sim);
	
	std::cout << "INFO: Initializing ILP ..." << std::endl;
	ClusteringILP2 ilp(sim, similarity_threshold_);
	ilp.initializeModel();
	
	std::cout << "INFO: Extending ILP ..." << std::endl;
	bool continue_computation = false;
     	for(size_t i=0; i<edges.size(); i+=1) {
		std::cout << "i: " << i+1 << " - " << std::get<2>(edges[i]) << std::endl;
		// Extend and check if cycles are introduced
		if(ilp.extendModel(edges, i, i+1)) continue;
		// Solve model
		continue_computation = ilp.solveModel(i+1);
		if(!continue_computation) break;

		// Save clustering
		auto clusters = TransitivityClusteringILP::ClusterEvaluation::extract_clusters(ilp.getSolution());
	        std::ofstream out;
        	out.open("out_" +  std::to_string(i+1) + ".txt");
	        for(size_t i=0; i<clusters.size(); ++i) {
        	        for(size_t j=0; j<clusters[i].size(); ++j) {
                	        out << sim.colName(clusters[i][j]) << "\t" << i << std::endl;
                	}
        	}
        	out.close();
	}

	return 0;
}

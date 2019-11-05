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
#include "IterativeClusteringILP.h"
#include "ClusteringILP.h"
#include "ClusterEvaluation.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string matrix_ = "", similarity_measure_ = "", output_prefix_ = "", output_edges_ = "", output_ = "";
bool isDistanceMatrix, runIteratively_, fixPreviousSolution_;
double similarity_threshold_, similarity_family_weight_;
size_t maxNumberOfIterations_, skipKIterations_, cplex_threads_, cplex_time_limit_;

MatrixReaderOptions matrixOptions;

bool parseArguments(int argc, char* argv[])
{
    bpo::variables_map vm;
    bpo::options_description desc;

    desc.add_options()("help,h", "Display this message")
        ("matrix,m", bpo::value<std::string>(&matrix_)->required(), "Name of a text file containing expression values as a matrix.")
        ("no-row-names,r", bpo::value<bool>(&matrixOptions.no_rownames)->default_value(false)->zero_tokens(), "The file does not contain row names.")
        ("no-col-names,c", bpo::value<bool>(&matrixOptions.no_colnames)->default_value(false)->zero_tokens(), "The file does not contain column names.")
        ("add-col-name,a", bpo::value<bool>(&matrixOptions.additional_colname)->default_value(false)->zero_tokens(), "File contains two lines specifying which rownames belong to which group.")
	
	("cplex-threads,t", bpo::value<size_t>(&cplex_threads_)->default_value(32), "Max number of threads that should be used by CPLEX.")
        ("cplex-time-limit,l", bpo::value<size_t>(&cplex_time_limit_)->default_value(10), "Time limit (in min) for each CPLEX iteration.")

        
        ("dissimilarity,d", bpo::value<bool>(&isDistanceMatrix)->default_value(false)->zero_tokens(), "The matrix already is a dissimilarity matrix.")
        ("similarity-measure,s", bpo::value<std::string>(&similarity_measure_)->default_value("pearson-correlation"), "Method used to compute pairwise dissimilarities.")
	
        ("run-itratively,i", bpo::value<bool>(&runIteratively_)->default_value(false)->zero_tokens(), "Run algorithm iteratively util threshold, maximum number of iteration, or time limit is reached.")
        ("fix-previous-solution,f", bpo::value<bool>(&fixPreviousSolution_)->default_value(false)->zero_tokens(), "Fix the solution of the previous iteration.")
        ("similarity-threshold,t", bpo::value<double>(&similarity_threshold_)->default_value(0.95), "Threshold for minimum dissimilarity score that should be used.")
	("max-number-of-iterations,x", bpo::value<size_t>(&maxNumberOfIterations_)->default_value(1000), "Max number of iterations.")
        ("skip-k-iterations,k", bpo::value<size_t>(&skipKIterations_)->default_value(0), "Number of iterations that should be performed before first model should be solved.")
        
        ("output-prefix,p", bpo::value<std::string>(&output_prefix_)->required(), "Prefix for all intermediate output files.")
        ("output-edges,e", bpo::value<std::string>(&output_edges_), "Output file for sorted edge.")
	("output,o", bpo::value<std::string>(&output_)->required(), "Name of the final output file.");
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
    } else if (similatity_measure == "angle-distance") {
        return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::AngleDistance()));
    } else if (similatity_measure == "normalized-angle-distance") {
        return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::NormalizedAngleDistance()));
    } else if (similatity_measure == "signed-distance-correlation") {
        return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::SignedDistanceCorrelation()));
    } else if (similatity_measure == "euclidean-distance") {
        return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::EuclideanDistance()));
    } else if (similatity_measure == "shifted-euclidean-distance-for-points") {
        return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::ShiftedEuclideanDistanceForPoints()));
    } else if (similatity_measure == "shifted-euclidean-distance-for-gradients") {
        return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::ShiftedEuclideanDistanceForGradients()));
    } else if (similatity_measure == "shifted-euclidean-distance-for-angles") {
        return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::ShiftedEuclideanDistanceForAngles()));
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
    if (output_edges_ != "") {
        std::ofstream out;
        out.open(output_edges_);
        for(auto edge : edges) {
            out << matrix.rowName(std::get<0>(edge)) << "\t" << matrix.rowName(std::get<1>(edge)) << "\t" << std::get<2>(edge) << std::endl; 
        }
        out.close();
    }
    return std::move(edges);
}

template<typename Matrix>
void saveSolution(IterativeClusteringILP& ilp, Matrix& sim, const std::string& file) {
    auto clusters = TransitivityClusteringILP::ClusterEvaluation::extract_clusters(ilp.getSolution());
    std::ofstream out;
    out.open(file);
    out << "class" << std::endl;
    for(size_t i=0; i<clusters.size(); ++i) {
        for(size_t j=0; j<clusters[i].size(); ++j) {
            out << sim.colName(clusters[i][j]) << "\t" << i+1 << std::endl;
        }
    }
    out.close();
}

template<typename Matrix>
void saveSolution(ClusteringILP& ilp, Matrix& sim, const std::string& file) {
    auto clusters = TransitivityClusteringILP::ClusterEvaluation::extract_clusters(ilp.getSolution());
    std::ofstream out;
    out.open(file);
    for(size_t i=0; i<clusters.size(); ++i) {
        for(size_t j=0; j<clusters[i].size(); ++j) {
            out << sim.colName(clusters[i][j]) << "\t" << i << std::endl;
        }
    }
    out.close();
}

template<typename Matrix>
void runIterativeClusteringILP(Matrix& sim) {
    std::cout << "INFO: Extracting sorted edges ..." << std::endl;
    std::vector<std::tuple<size_t, size_t, double>> edges = extract_sorted_edges(sim);
        
    std::cout << "INFO: Initializing ILP ..." << std::endl;
    IterativeClusteringILP ilp(sim, cplex_threads_, cplex_time_limit_);
    ilp.initializeModel();
       
    std::cout << "INFO: Extending ILP ..." << std::endl;
    bool continue_computation = false;
    for(size_t i=0; i<edges.size(); i+=1) {
        std::cout << "i: " << i+1 << " - " << std::get<2>(edges[i]) << std::endl;
	bool stop = similarity_threshold_ > std::get<2>(edges[i]) || i > maxNumberOfIterations_;
	// Extend model and continue if no constraints are introduced
        if(ilp.extendModel(edges, i, i+1, fixPreviousSolution_) || (i < skipKIterations_) || (!runIteratively_ && !stop)) continue;
        // Solve model and stop computation if model is infeasible or time limit was exceeded
        continue_computation = ilp.solveModel() && !stop;
        if(!continue_computation) break;
        // Save last feasible solution
        saveSolution(ilp , sim, output_prefix_ +  std::to_string(i+1) + ".txt");
    }
    saveSolution(ilp , sim, output_);
}

/*template<typename Matrix>
void runThresholdClusteringILP(Matrix& sim) {        
    std::cout << "INFO: Initializing ILP ..." << std::endl;
    ClusteringILP ilp(sim, similarity_threshold_, ClusteringILP::Metric::similarity);
    ilp.buildModel();
    saveSolution(ilp , sim, "out_t_" + std::to_string(similarity_threshold_) + ".txt");
}*/

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
    try {
        matrix = readDenseMatrix(matrix_, matrixOptions);
    } catch(const IOError& e) {
        std::cerr << "ERROR: Could not open input data matrix for reading." << std::endl;
        return -4;
    }

    if(!isDistanceMatrix) {
        std::cout << "INFO: Calculating similatiry matrix ..." << std::endl;
        DenseMatrix sim = compute_similarity_matrix(matrix, similarity_measure_);
        runIterativeClusteringILP(sim);
    } else {
        runIterativeClusteringILP(matrix);
    } 

    return 0;
}

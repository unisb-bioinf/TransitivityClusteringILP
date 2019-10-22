#include <iostream>

#include <boost/program_options.hpp>

#include "DenseMatrixReader.h"
#include "DenseMatrix.h"
#include "matrixTools.h"
#include "MatrixIterator.h"
#include "SimilarityMeasures.h"

using namespace GeneTrail;
namespace bpo = boost::program_options;

std::string matrix_ = "", similarity_measure_ = "";

MatrixReaderOptions matrixOptions;

bool parseArguments(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	desc.add_options()("help,h", "Display this message")
		("matrix,m", bpo::value<std::string>(&matrix_)->required(), "Name of a text file containing expression values as a matrix.")
		("no-row-names,r", bpo::value<bool>(&matrixOptions.no_rownames)->default_value(false)->zero_tokens(), "Does the file contain row names.")
		("no-col-names,c", bpo::value<bool>(&matrixOptions.no_colnames)->default_value(false)->zero_tokens(), "Does the file contain column names.")
		("add-col-name,a", bpo::value<bool>(&matrixOptions.additional_colname)->default_value(false)->zero_tokens(), "File containing two lines specifying which rownames belong to which group.")
		("similarity-measure,s", bpo::value<std::string>(&similarity_measure_)->required(), "Method used to compute pairwise similarities.");
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
	for(size_t i=0; i<matrix.rows(); ++i){
		dist.set(i, i, 1.0);
		for(size_t j=i+1; i<matrix.rows(); ++i){
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
	} else if(similatity_measure == "spearman-correlation") {
                return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::SpearmanCorrelation()));
        } else {
                return std::move(compute_similarity_matrix(matrix, TransitivityClusteringILP::ShiftedEuclideanDistanceForGradientsAndPoints()));
        }
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

	try {
		matrix = readDenseMatrix(matrix_, matrixOptions);
	} catch(const IOError& e) {
		std::cerr << "ERROR: Could not open input data matrix for reading." << std::endl;
		return -4;
	}

	std::cout << "INFO: Calculating similatiry matrix" << std::endl;
	compute_similarity_matrix(matrix, similarity_measure_);

	return 0;
}

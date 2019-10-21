#ifndef MATRIX_TOOLS_H
#define MATRIX_TOOLS_H

#include "DenseColumnSubset.h"

#include "macros.h"

#include <string>
#include <vector>

namespace GeneTrail
{
	class DenseMatrix;

	struct GT2_EXPORT MatrixReaderOptions
	{
		bool no_rownames = false;
		bool no_colnames = false;
		bool additional_colname = false;
	};

	class GT2_EXPORT EmptyGroup : public std::exception
	{
		public:
		EmptyGroup(const std::string& name) noexcept : groupname_(name) {}

		virtual const char* what() const noexcept
		{
			return (std::string("Group \"") + groupname_ +
			        "\" does not contain any datapoint.").c_str();
		}

		private:
		std::string groupname_;
	};

	GT2_EXPORT std::tuple<DenseColumnSubset, DenseColumnSubset>
	splitMatrix(DenseMatrix& matrix, const std::vector<std::string>& reference,
	            const std::vector<std::string>& test);
	GT2_EXPORT DenseMatrix buildDenseMatrix(const std::string& expr1,
	                                        const std::string& expr2,
	                                        const MatrixReaderOptions& options);
	GT2_EXPORT DenseMatrix readDenseMatrix(const std::string& matrix,
	                                       const MatrixReaderOptions& options);
	std::vector<unsigned int> getIndices(const DenseMatrix& matrix, const std::vector<std::string>& colnames, const std::string& groupname);
}

#endif // MATRIX_TOOLS_H

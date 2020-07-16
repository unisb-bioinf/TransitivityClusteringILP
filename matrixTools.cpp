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

#include "matrixTools.h"

#include "DenseMatrix.h"
#include "DenseMatrixReader.h"

#include <iostream>
#include <fstream>

namespace GeneTrail
{

	std::vector<unsigned int>
	getIndices(const DenseMatrix& matrix,
	           const std::vector<std::string>& colnames,
	           const std::string& groupname)
	{
		std::vector<unsigned int> indices;
		indices.reserve(colnames.size());
		for(const auto& s : colnames) {
			if(matrix.hasCol(s)) {
				indices.emplace_back(matrix.colIndex(s));
			} else {
				std::cerr << "WARNING: Could not find column \"" + s + "\".\n";
			}
		}

		if(indices.empty() && !colnames.empty()) {
			throw EmptyGroup(groupname);
		}

		return indices;
	}

	std::tuple<DenseColumnSubset, DenseColumnSubset>
	splitMatrix(DenseMatrix& matrix, const std::vector<std::string>& reference,
	            const std::vector<std::string>& test)
	{
		return std::make_tuple(
		    DenseColumnSubset(
		        &matrix, getIndices(matrix, reference, "reference")),
		    DenseColumnSubset(
		        &matrix, getIndices(matrix, test, "test")));
	}

	DenseMatrix buildDenseMatrix(const std::string& expr1,
	                             const std::string& expr2,
	                             const MatrixReaderOptions& options)
	{
		auto m1 = readDenseMatrix(expr1, options);
		if(expr2 != "") {
			auto m2 = readDenseMatrix(expr2, options);
			m1.cbind(m2);
		}
		return m1;
	}

	DenseMatrix readDenseMatrix(const std::string& matrix,
	                            const MatrixReaderOptions& options)
	{
		unsigned int opts = DenseMatrixReader::NO_OPTIONS;

		if(!options.no_rownames) {
			opts |= DenseMatrixReader::READ_ROW_NAMES;
		}

		if(!options.no_colnames) {
			opts |= DenseMatrixReader::READ_COL_NAMES;
		}

		if(options.additional_colname) {
			opts |= DenseMatrixReader::ADDITIONAL_COL_NAME;
		}

		DenseMatrixReader reader;

		std::ifstream strm(matrix, std::ios::binary);
		return reader.read(strm, opts);
	}
}

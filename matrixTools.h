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

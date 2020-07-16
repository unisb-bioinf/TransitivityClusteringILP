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

#ifndef GT2_DENSE_MATRIX_WRITER_H
#define GT2_DENSE_MATRIX_WRITER_H

#include "macros.h"
#include "MatrixWriter.h"

#include <ostream>
#include <vector>

namespace GeneTrail
{
	class Matrix;
	class DenseMatrix;

	class GT2_EXPORT DenseMatrixWriter : public MatrixWriter
	{
		public:
			void writeText(std::ostream& output, const Matrix& matrix) const;

			/**
			 * Reads a matrix from a binary file.
			 * 
			 * \see DenseMatrixReader::binaryRead_
			 */
			uint64_t writeBinary(std::ostream& output, const DenseMatrix& matrix) const;
			uint64_t writeBinary(std::ostream& output, const Matrix& matrix) const;

		private:
			uint64_t writeData_(std::ostream& output, const DenseMatrix& matrix) const;
			uint64_t writeData_(std::ostream& output, const Matrix& matrix) const;
	};
}

#endif //GT2_DENSE_MATRIX_WRITER_H


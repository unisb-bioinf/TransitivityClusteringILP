/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2016 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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

#ifndef GT2_DENSE_ROW_SUBSET_H
#define GT2_DENSE_ROW_SUBSET_H

#include "DenseMatrix.h"

#include "macros.h"

namespace GeneTrail
{
	/**
	 * This class represents a subset of a DenseMatrix.
	 *
	 * The reason we derive from DenseMatrix is that we
	 * want to be able to use DenseMatrix and DenseRowSubset
	 * interchangably.
	 *
	 * Lets see if this works...
	 */
	class GT2_EXPORT DenseRowSubset : public Matrix
	{
		public:
			using DMatrix = DenseMatrix::DMatrix;
			typedef std::vector<DenseMatrix::index_type> ISubset;
			typedef std::vector<std::string> SSubset;

			template<typename Iterator>
			DenseRowSubset(DenseMatrix* mat, Iterator first, Iterator last);
			DenseRowSubset(DenseMatrix* mat, ISubset rows);
			DenseRowSubset(DenseMatrix* mat, const SSubset& rows);

			DenseRowSubset(const DenseRowSubset& subs) = default;
			DenseRowSubset(DenseRowSubset&& subs);

			DenseRowSubset& operator=(const DenseRowSubset& subs) = default;
			DenseRowSubset& operator=(DenseRowSubset&& subs);

			template<typename InputIterator>
			void assign(InputIterator first, InputIterator last);

			value_type& operator()(index_type i, index_type j) override
			{
				return (*mat_)(row_subset_[i], j);
			}

			value_type operator()(index_type i, index_type j) const override
			{
				return (*mat_)(row_subset_[i], j);
			}

			/**
			 * Direct access to the Eigen expression template representing the
			 * respective row. This allows leveraging Eigen's vectorization
			 * features.
			 */
			DMatrix::RowXpr row(index_type i) {
				return mat_->matrix().row(row_subset_[i]);
			}

			/**
			 * Direct access to the Eigen expression template representing the
			 * respective row. This allows leveraging Eigen's vectorization
			 * features.
			 */
			DMatrix::ConstRowXpr row(index_type i) const {
				return static_cast<const DenseMatrix*>(mat_)->matrix().row(row_subset_[i]);
			}

			const std::string& colName(index_type j) const override;
			const std::string& rowName(index_type i) const override;

			index_type colIndex(const std::string& col) const override;
			index_type rowIndex(const std::string& row) const override;

			index_type cols() const override;
			index_type rows() const override;

			bool hasCol(const std::string& name) const override;
			bool hasRow(const std::string& name) const override;

			void setColName(index_type j, const std::string& new_name) override;
			void setColName(const std::string& old_name, const std::string& new_name) override;
			void setColNames(const std::vector< std::string >& col_names) override;
			void setRowName(index_type i, const std::string& new_name) override;
			void setRowName(const std::string& old_name, const std::string& new_name) override;
			void setRowNames(const std::vector< std::string >& row_names) override;

			void removeCols(const std::vector< index_type >& indices) override;
			void removeRows(const std::vector< index_type >& indices) override;
			void shuffleCols(const std::vector< index_type >& perm) override;
			void shuffleRows(const std::vector< index_type >& perm) override;
			void transpose() override;

			const std::vector< std::string >& colNames() const override;
			const std::vector< std::string >& rowNames() const override;

		private:
			DenseMatrix* mat_;
			ISubset row_subset_;

			// TODO: Think of something smart to fix the hack below
			mutable SSubset row_names_cache_;

			void remove_(const std::vector<Matrix::index_type>& indices, ISubset& subset);
	};

	template<typename Iterator>
	DenseRowSubset::DenseRowSubset(DenseMatrix* mat, Iterator begin, Iterator end)
		: mat_(mat),
		  row_subset_(begin, end)
	{
	}

	template<typename InputIterator>
	void DenseRowSubset::assign(InputIterator first, InputIterator last)
	{
		row_subset_.assign(first, last);
	}

}

#endif // GT2_DENSE_MATRIX_SUBSET_H


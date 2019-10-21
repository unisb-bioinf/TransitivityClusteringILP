/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2013-2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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

#include "DenseRowSubset.h"

#include "Exception.h"

#include <algorithm>

namespace GeneTrail
{
	/*
	 * Constructors
	 */
	DenseRowSubset::DenseRowSubset(DenseMatrix* mat, const SSubset& rows)
		: mat_(mat),
		  row_subset_(rows.size())
	{
		std::transform(rows.begin(), rows.end(), row_subset_.begin(), [mat](const std::string& s) {
			index_type i = mat->rowIndex(s);

			if(i == std::numeric_limits<index_type>::max()) {
				throw InvalidKey(s);
			}

			return i;
		});
	}

	DenseRowSubset::DenseRowSubset(DenseMatrix* mat, ISubset rows)
		: mat_(mat),
		  row_subset_(std::move(rows))
	{
	}

	DenseRowSubset::DenseRowSubset(DenseRowSubset&& subs)
		: mat_(subs.mat_),
		  row_subset_(std::move(subs.row_subset_))
	{
	}

	/*
	 * Operators
	 */
	DenseRowSubset& DenseRowSubset::operator=(DenseRowSubset&& subs)
	{
		assert(this != &subs);

		mat_ = subs.mat_;
		row_subset_ = std::move(subs.row_subset_);

		return *this;
	}

	Matrix::index_type DenseRowSubset::cols() const
	{
		return mat_->cols();
	}

	Matrix::index_type DenseRowSubset::rows() const
	{
		return row_subset_.size();
	}

	bool DenseRowSubset::hasCol(const std::string& name) const
	{
		return mat_->hasCol(name);
	}

	bool DenseRowSubset::hasRow(const std::string& name) const
	{
		for(auto i : row_subset_) {
			if(mat_->rowName(i) == name) {
				return true;
			}
		}

		return false;
	}

	Matrix::index_type DenseRowSubset::rowIndex(const std::string& row) const
	{
		int idx = 0;

		for(auto i : row_subset_) {

			if(mat_->rowName(i) == row) {
				return idx;
			}
			++idx;
		}

		return -1;
	}

	const std::string& DenseRowSubset::rowName(Matrix::index_type i) const
	{
		return mat_->rowName(row_subset_[i]);
	}

	const std::vector< std::string >& DenseRowSubset::rowNames() const
	{
		row_names_cache_.resize(row_subset_.size());

		for(size_t i = 0; i < row_subset_.size(); ++i) {
			row_names_cache_[i] = mat_->rowName(row_subset_[i]);
		}

		return row_names_cache_;
	}

	Matrix::index_type DenseRowSubset::colIndex(const std::string& col) const
	{
		return mat_->colIndex(col);
	}

	const std::string& DenseRowSubset::colName(Matrix::index_type j) const
	{
		return mat_->colName(j);
	}

	const std::vector< std::string >& DenseRowSubset::colNames() const
	{
		return mat_->colNames();
	}

	void DenseRowSubset::setColName(Matrix::index_type j, const std::string& new_name)
	{
		mat_->setColName(j, new_name);
	}

	void DenseRowSubset::setColName(const std::string& old_name, const std::string& new_name)
	{
		setColName(colIndex(old_name), new_name);
	}

	void DenseRowSubset::setColNames(const std::vector< std::string >& col_names)
	{
		int i = 0;
		for(const auto& s : col_names) {
			mat_->setColName(i, s);
			++i;
		}
	}

	void DenseRowSubset::setRowName(Matrix::index_type i, const std::string& new_name)
	{
		mat_->setRowName(row_subset_[i], new_name);
	}

	void DenseRowSubset::setRowName(const std::string& old_name, const std::string& new_name)
	{
		setRowName(rowIndex(old_name), new_name);
	}

	void DenseRowSubset::setRowNames(const std::vector< std::string >& row_names)
	{
		int i = 0;
		for(const auto& s : row_names) {
			mat_->setRowName(row_subset_[i], s);
			++i;
		}
	}

	/*
	 * Matrix Operations
	 */
	void DenseRowSubset::remove_(const std::vector< Matrix::index_type >& indices, ISubset& subset)
	{
		size_t next_idx = 1;
		size_t write_idx = indices[0];

		for(size_t read_idx = indices[0] + 1; read_idx < subset.size(); ++read_idx) {
			if(next_idx < indices.size() && read_idx == indices[next_idx]) {
				assert(indices[next_idx -1] < indices[next_idx]);
				++next_idx;
			}
			else
			{
				subset[write_idx] = subset[read_idx];
				++write_idx;
			}
		}

		subset.resize(subset.size() - indices.size());
	}

	void DenseRowSubset::removeCols(const std::vector< Matrix::index_type >& indices)
	{
		mat_->removeCols(indices);
	}

	void DenseRowSubset::removeRows(const std::vector< Matrix::index_type >& indices)
	{
		remove_(indices, row_subset_);
	}

	void DenseRowSubset::shuffleCols(const std::vector< Matrix::index_type >& perm)
	{
		mat_->shuffleCols(perm);
	}

	void DenseRowSubset::shuffleRows(const std::vector< Matrix::index_type >& perm)
	{
		ISubset tmp(row_subset_.size());

		for(size_t i = 0; i < perm.size(); ++i) {
			tmp[i] = row_subset_[perm[i]];
		}

		std::swap(row_subset_, tmp);
	}

	void DenseRowSubset::transpose()
	{
		throw NotImplemented(__FILE__, __LINE__, "DenseRowSubset::transpose()");
	}

}


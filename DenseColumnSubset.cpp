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

#include "DenseColumnSubset.h"

#include "Exception.h"

#include <algorithm>

namespace GeneTrail
{
	/*
	 * Constructors
	 */
	DenseColumnSubset::DenseColumnSubset(DenseMatrix* mat, const SSubset& cols)
		: mat_(mat),
		  col_subset_(cols.size())
	{
		std::transform(cols.begin(), cols.end(), col_subset_.begin(), [mat](const std::string& s) {
			index_type i = mat->colIndex(s);

			if(i == std::numeric_limits<index_type>::max()) {
				throw InvalidKey(s);
			}

			return i;
		});
	}

	DenseColumnSubset::DenseColumnSubset(DenseMatrix* mat, ISubset cols)
		: mat_(mat),
		  col_subset_(std::move(cols))
	{
	}

	DenseColumnSubset::DenseColumnSubset(DenseColumnSubset&& subs)
		: mat_(subs.mat_),
		  col_subset_(std::move(subs.col_subset_))
	{
	}

	/*
	 * Operators
	 */
	DenseColumnSubset& DenseColumnSubset::operator=(DenseColumnSubset&& subs)
	{
		assert(this != &subs);

		mat_ = subs.mat_;
		col_subset_ = std::move(subs.col_subset_);

		return *this;
	}

	Matrix::index_type DenseColumnSubset::cols() const
	{
		return col_subset_.size();
	}

	Matrix::index_type DenseColumnSubset::rows() const
	{
		return mat_->rows();
	}

	bool DenseColumnSubset::hasRow(const std::string& name) const
	{
		return mat_->hasRow(name);
	}

	bool DenseColumnSubset::hasCol(const std::string& name) const
	{
		for(auto i : col_subset_) {
			if(mat_->colName(i) == name) {
				return true;
			}
		}

		return false;
	}

	Matrix::index_type DenseColumnSubset::rowIndex(const std::string& row) const
	{
		return mat_->rowIndex(row);
	}

	const std::string& DenseColumnSubset::rowName(Matrix::index_type i) const
	{
		return mat_->rowName(i);
	}

	const std::vector< std::string >& DenseColumnSubset::rowNames() const
	{
		return mat_->rowNames();
	}

	Matrix::index_type DenseColumnSubset::colIndex(const std::string& col) const
	{
		int idx = 0;

		for(auto i : col_subset_) {

			if(mat_->colName(i) == col) {
				return idx;
			}
			++idx;
		}

		return -1;
	}

	const std::string& DenseColumnSubset::colName(Matrix::index_type j) const
	{
		return mat_->colName(col_subset_[j]);
	}

	const std::vector< std::string >& DenseColumnSubset::colNames() const
	{
		col_names_cache_.resize(col_subset_.size());

		for(size_t i = 0; i < col_subset_.size(); ++i) {
			col_names_cache_[i] = mat_->colName(col_subset_[i]);
		}

		return col_names_cache_;
	}

	void DenseColumnSubset::setColName(Matrix::index_type j, const std::string& new_name)
	{
		mat_->setColName(col_subset_[j], new_name);
	}

	void DenseColumnSubset::setColName(const std::string& old_name, const std::string& new_name)
	{
		setColName(colIndex(old_name), new_name);
	}

	void DenseColumnSubset::setColNames(const std::vector< std::string >& col_names)
	{
		int i = 0;
		for(const auto& s : col_names) {
			mat_->setColName(col_subset_[i], s);
			++i;
		}
	}

	void DenseColumnSubset::setRowName(Matrix::index_type i, const std::string& new_name)
	{
		mat_->setRowName(i, new_name);
	}

	void DenseColumnSubset::setRowName(const std::string& old_name, const std::string& new_name)
	{
		setRowName(rowIndex(old_name), new_name);
	}

	void DenseColumnSubset::setRowNames(const std::vector< std::string >& row_names)
	{
		mat_->setRowNames(row_names);
	}

	/*
	 * Matrix Operations
	 */
	void DenseColumnSubset::remove_(const std::vector< Matrix::index_type >& indices, ISubset& subset)
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

	void DenseColumnSubset::removeCols(const std::vector< Matrix::index_type >& indices)
	{
		remove_(indices, col_subset_);
	}

	void DenseColumnSubset::removeRows(const std::vector< Matrix::index_type >& indices)
	{
		mat_->removeRows(indices);
	}

	void DenseColumnSubset::shuffleCols(const std::vector< Matrix::index_type >& perm)
	{
		ISubset tmp(col_subset_.size());

		for(size_t i = 0; i < perm.size(); ++i) {
			tmp[i] = col_subset_[perm[i]];
		}

		std::swap(col_subset_, tmp);
	}

	void DenseColumnSubset::shuffleRows(const std::vector< Matrix::index_type >& perm)
	{
		mat_->shuffleRows(perm);
	}

	void DenseColumnSubset::transpose()
	{
		throw NotImplemented(__FILE__, __LINE__, "DenseColumnSubset::transpose()");
	}

}


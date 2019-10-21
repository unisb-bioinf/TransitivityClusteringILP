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

#ifndef GT2_CORE_DENSEMATRIXITERATOR_H
#define GT2_CORE_DENSEMATRIXITERATOR_H

#include <cassert>
#include <cstddef>
#include <iterator>

namespace GeneTrail
{

template <typename Matrix, typename T, typename Access>
class DenseMatrixIterator
    : std::iterator<std::random_access_iterator_tag, T>
{
  public:
	using value_type = T;

	explicit DenseMatrixIterator(const Matrix* matrix = nullptr,
	                           size_t index = 0) noexcept
		: matrix_(matrix),
		  index_(index)
	{
	}

	DenseMatrixIterator(const DenseMatrixIterator&) noexcept = default;
	DenseMatrixIterator(DenseMatrixIterator&&) noexcept = default;

	DenseMatrixIterator& operator=(const DenseMatrixIterator&) noexcept = default;
	DenseMatrixIterator& operator=(DenseMatrixIterator&&) noexcept = default;

	bool operator==(const DenseMatrixIterator& o) const noexcept
	{
		return matrix_ == o.matrix_ && index_ == o.index_;
	}

	bool operator!=(const DenseMatrixIterator& o) const noexcept
	{
		return matrix_ != o.matrix_ || index_ != o.index_;
	}

	value_type operator*() const { return Access()(matrix_, index_); }

	value_type operator->() { return Access()(matrix_, index_); }

	const value_type& operator->() const { return Access()(matrix_, index_); }

	value_type operator[](size_t n) { return Access()(matrix_, index_ + n); }

	DenseMatrixIterator& operator++()
	{
		++index_;
		return *this;
	}

	DenseMatrixIterator operator++(int)
	{
		auto tmp = index_++;
		return DenseMatrixIterator(matrix_, tmp);
	}

	DenseMatrixIterator& operator--()
	{
		--index_;
		return *this;
	}

	DenseMatrixIterator operator--(int)
	{
		auto tmp = index_--;
		return DenseMatrixIterator(matrix_, tmp);
	}

	DenseMatrixIterator operator+(std::ptrdiff_t n) const
	{
		return DenseMatrixIterator(matrix_, index_ + n);
	}

	DenseMatrixIterator operator-(std::ptrdiff_t n) const
	{
		return DenseMatrixIterator(matrix_, index_ - n);
	}

	std::ptrdiff_t operator-(const DenseMatrixIterator& o) const
	{
		return index_ - o.index_;
	}

	bool operator<(const DenseMatrixIterator& o) const
	{
		assert(matrix_ == o.matrix_);
		return index_ < o.index_;
	}

	bool operator<=(const DenseMatrixIterator& o) const
	{
		assert(matrix_ == o.matrix_);
		return index_ <= o.index_;
	}

	bool operator>(const DenseMatrixIterator& o) const
	{
		assert(matrix_ == o.matrix_);
		return index_ > o.index_;
	}

	bool operator>=(const DenseMatrixIterator& o) const
	{
		assert(matrix_ == o.matrix_);
		return index_ >= o.index_;
	}

	void operator+=(std::ptrdiff_t n) { index_ += n; }

	void operator-=(std::ptrdiff_t n) { index_ -= n; }

  private:
	const Matrix* matrix_;
	size_t index_;
};

template <typename M, typename T, typename A>
DenseMatrixIterator<M, T, A> operator+(std::ptrdiff_t n, const DenseMatrixIterator<M, T, A>& it)
{
	return it + n;
}

template<typename Matrix>
struct ColumnAccess {
	typename Matrix::DMatrix::ConstColXpr operator()(const Matrix* m, size_t index) const {
		return m->col(index);
	}
};

template<typename Matrix>
struct RowAccess {
	typename Matrix::DMatrix::ConstRowXpr operator()(const Matrix* m, size_t index) const {
		return m->row(index);
	}
};

template<typename Matrix>
using DenseColumnIterator = DenseMatrixIterator<Matrix, typename Matrix::DMatrix::ConstColXpr, ColumnAccess<Matrix>>;

template<typename Matrix>
using DenseRowIterator = DenseMatrixIterator<Matrix, typename Matrix::DMatrix::ConstRowXpr, RowAccess<Matrix>>;

template<typename Matrix>
DenseColumnIterator<Matrix> make_column_iterator(const Matrix* m, size_t i)
{
	return DenseColumnIterator<Matrix>(m, i);
}

template<typename Matrix>
DenseRowIterator<Matrix> make_row_iterator(const Matrix* m, size_t i)
{
	return DenseRowIterator<Matrix>(m, i);
}

}

#endif // GT2_CORE_DENSEMATRIXITERATOR_H

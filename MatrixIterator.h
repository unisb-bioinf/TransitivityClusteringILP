/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
#ifndef GT2_MATRIX_ITERATOR_H
#define GT2_MATRIX_ITERATOR_H

#include <cstddef>
#include <iterator>

namespace GeneTrail
{
	template <typename Matrix>
	class ColumnIterator : public std::iterator<std::random_access_iterator_tag,
	                                            typename Matrix::value_type,
	                                            std::ptrdiff_t,
	                                            typename Matrix::value_type*,
	                                            typename Matrix::value_type>
	{
		public:
		using index_type = typename Matrix::index_type;
		using value_type = typename Matrix::value_type;

		ColumnIterator(const Matrix* matrix, index_type row, index_type col)
		    : matrix_(matrix), row_(row), col_(col)
		{
		}

		ColumnIterator(const ColumnIterator&) = default;
		ColumnIterator& operator=(const ColumnIterator&) = default;

		bool operator==(const ColumnIterator& it) const
		{
			return col_ == it.col_ && row_ == it.row_ && matrix_ == it.matrix_;
		}

		bool operator!=(const ColumnIterator& it) const
		{
			return col_ != it.col_ || row_ != it.row_ || matrix_ != it.matrix_;
		}

		value_type operator*() const { return (*matrix_)(row_, col_); }

		ColumnIterator& operator++()
		{
			++col_;
			return *this;
		}

		ColumnIterator operator++(int)
		{
			ColumnIterator tmp(*this);
			++col_;
			return tmp;
		}

		ColumnIterator& operator--()
		{
			--col_;
			return *this;
		}

		ColumnIterator operator--(int)
		{
			ColumnIterator tmp(*this);
			--col_;
			return tmp;
		}

		template <typename T>
		friend ColumnIterator<T> operator+(const ColumnIterator<T>& a,
		                                   std::ptrdiff_t n);
		template <typename T>
		friend ColumnIterator<T> operator+(std::ptrdiff_t n,
		                                   const ColumnIterator<T>& a);
		template <typename T>
		friend ColumnIterator<T> operator-(const ColumnIterator<T>& a,
		                                   std::ptrdiff_t n);
		template <typename T>
		friend std::ptrdiff_t operator-(const ColumnIterator<T>& a,
		                                const ColumnIterator<T>& b);

		bool operator<(const ColumnIterator& it) const
		{
			return col_ < it.col_;
		}

		bool operator>(const ColumnIterator& it) const
		{
			return col_ > it.col_;
		}

		bool operator<=(const ColumnIterator& it) const
		{
			return col_ <= it.col_;
		}

		bool operator>=(const ColumnIterator& it) const
		{
			return col_ >= it.col_;
		}

		void operator+=(std::ptrdiff_t n) { col_ += n; }

		void operator-=(std::ptrdiff_t n) { col_ -= n; }

		value_type operator[](std::ptrdiff_t n) const
		{
			return (*matrix_)(row_, col_ + n);
		}

		index_type row() const { return row_; }
		index_type col() const { return col_; }

		private:
		const Matrix* matrix_;
		index_type row_;
		index_type col_;
	};

	template <typename Matrix>
	ColumnIterator<Matrix> operator+(const ColumnIterator<Matrix>& a,
	                                 std::ptrdiff_t n)
	{
		return ColumnIterator<Matrix>(a.matrix_, a.row_, a.col_ + n);
	}

	template <typename Matrix>
	ColumnIterator<Matrix> operator+(std::ptrdiff_t n,
	                                 ColumnIterator<Matrix>& a)
	{
		return ColumnIterator<Matrix>(a.matrix_, a.row_, a.col_ + n);
	}

	template <typename Matrix>
	ColumnIterator<Matrix> operator-(const ColumnIterator<Matrix>& a,
	                                 std::ptrdiff_t n)
	{
		return ColumnIterator<Matrix>(a.matrix_, a.row_, a.col_ - n);
	}

	template <typename Matrix>
	std::ptrdiff_t operator-(const ColumnIterator<Matrix>& a,
	                         const ColumnIterator<Matrix>& b)
	{
		return a.col_ - b.col_;
	}

	template <typename Matrix> class ColumnProxy
	{
		using index_type = typename Matrix::index_type;

		public:
		ColumnProxy(const Matrix* matrix, index_type row)
		    : matrix_(matrix), row_(row)
		{
		}
		ColumnIterator<Matrix> begin() const
		{
			return ColumnIterator<Matrix>(matrix_, row_, 0);
		}

		ColumnIterator<Matrix> end() const
		{
			return ColumnIterator<Matrix>(matrix_, row_, matrix_->cols());
		}

		private:
		const Matrix* matrix_;
		index_type row_;
	};

	template <typename Matrix>
	class RowMajorMatrixIterator
	    : public std::iterator<std::forward_iterator_tag, ColumnProxy<Matrix>>
	{
		public:
		using index_type = typename Matrix::index_type;
		using value_type = ColumnProxy<Matrix>;

		RowMajorMatrixIterator(const Matrix* matrix, index_type row)
		    : proxy_(matrix, row), matrix_(matrix), row_(row)
		{
		}

		RowMajorMatrixIterator(const RowMajorMatrixIterator&) = default;
		RowMajorMatrixIterator&
		operator=(const RowMajorMatrixIterator&) = default;

		bool operator==(const RowMajorMatrixIterator& it)
		{
			return matrix_ == it.matrix_ && row_ == it.row_;
		}

		bool operator!=(const RowMajorMatrixIterator& it)
		{
			return !(*this == it);
		}

		value_type operator*() const { return value_type(matrix_, row_); }

		value_type* operator->() const
		{
			proxy_ = ColumnProxy<Matrix>(matrix_, row_);

			return &proxy_;
		}

		RowMajorMatrixIterator& operator++()
		{
			++row_;
			return *this;
		}

		RowMajorMatrixIterator operator++(int)
		{
			RowMajorMatrixIterator tmp(*this);
			++row_;
			return tmp;
		}

		index_type row() const { return row_; }

		private:
		mutable ColumnProxy<Matrix> proxy_;
		const Matrix* matrix_;
		index_type row_;
	};
}

#endif // GT2_MATRIX_ITERATOR_H

/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
 * for more details.                                                            *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/

#ifndef CGOGN_MODELING_ALGOS_DECIMATION_CELL_QUEUE_H_
#define CGOGN_MODELING_ALGOS_DECIMATION_CELL_QUEUE_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/types/mesh_traits.h>

#include <map>

namespace cgogn
{

namespace modeling
{

template <typename CELL>
class CellQueue
{
public:
	using Self = CellQueue<CELL>;
	using EdgeCrit = std::multimap<cgogn::float64, CELL>;

	struct CellQueueInfo
	{
		typename EdgeCrit::const_iterator it_;
		bool valid_;
		CellQueueInfo() : valid_(false)
		{
		}
	};

	inline CellQueue()
	{
	}
	virtual ~CellQueue()
	{
	}

	// set to 64 bit to avoid conversion warning, can be 32 but need cast on insertion
	std::multimap<cgogn::float64, CELL> cells_;

	class const_iterator
	{
	public:
		const Self* const queue_ptr_;
		typename EdgeCrit::const_iterator cell_it_;

		inline const_iterator(const Self* trav, typename EdgeCrit::const_iterator it) : queue_ptr_(trav), cell_it_(it)
		{
		}
		inline const_iterator(const const_iterator& it) : queue_ptr_(it.queue_ptr_), cell_it_(it.cell_it_)
		{
		}

		inline const_iterator& operator=(const const_iterator& it)
		{
			queue_ptr_ = it.queue_ptr_;
			cell_it_ = it.cell_it_;
			return *this;
		}
		inline const_iterator& operator++()
		{
			cell_it_ = queue_ptr_->cells_.begin();
			return *this;
		}
		inline const CELL& operator*() const
		{
			return (*cell_it_).second;
		}
		inline bool operator!=(const_iterator it) const
		{
			cgogn_assert(queue_ptr_ == it.queue_ptr_);
			return cell_it_ != it.cell_it_;
		}
	};

	const_iterator begin() const
	{
		return const_iterator(this, cells_.begin());
	}
	const_iterator end() const
	{
		return const_iterator(this, cells_.end());
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_DECIMATION_CELL_QUEUE_H_

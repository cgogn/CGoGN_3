/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_CORE_CONTAINER_VECTOR_H_
#define CGOGN_CORE_CONTAINER_VECTOR_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/numerics.h>

#include <cgogn/core/types/container/attribute_container.h>

#include <vector>
#include <string>
#include <memory>

namespace cgogn
{

//////////////////
// Vector class //
//////////////////

template <typename T>
class CGOGN_CORE_EXPORT Vector : public AttributeGen
{
private:

	std::vector<T> data_;

	void manage_index(uint32 index) override
	{
		while (index >= data_.size())
			data_.push_back(T());
	}

public:

	Vector(AttributeContainerGen* container, bool is_mark, const std::string& name) : AttributeGen(container, is_mark, name)
	{
		data_.reserve(512u);
	}

	~Vector() override
	{}

	inline T& operator[](uint32 index)
	{
		cgogn_message_assert(index < data_.size(), "index out of bounds");
		return data_[index];
	}

	inline const T& operator[](uint32 index) const
	{
		cgogn_message_assert(index < data_.size(), "index out of bounds");
		return data_[index];
	}

	inline void swap(Vector<T>* ca)
	{
		if (ca->container_ == this->container_)
			data_.swap(ca->data_);
	}

	inline const void* data_pointer() const
	{
		return &data_[0];
	}

	class const_iterator
	{
		const Vector<T>* ca_;
		uint32 index_;

	public:

		inline const_iterator(const Vector<T>* ca, uint32 index) : ca_(ca), index_(index)
		{}
		inline const_iterator(const const_iterator& it) : ca_(it.ca_), index_(it.index_)
		{}
		inline const_iterator& operator=(const const_iterator& it)
		{
			ca_ = it.ca_;
			index_ = it.index_;
			return *this;
		}
		inline bool operator!=(const_iterator it) const
		{
			cgogn_assert(ca_ == it.ca_);
			return index_ != it.index_;
		}
		inline const_iterator& operator++()
		{
			index_ = ca_->container_->next_index(index_);
			return *this;
		}
		inline const T& operator*() const
		{
			return ca_->operator[](index_);
		}
		inline uint32 index() { return index_; }
	};
	inline const_iterator begin() const { return const_iterator(this, this->container_->first_index()); }
	inline const_iterator end() const { return const_iterator(this, this->container_->last_index()); }

	class iterator
	{
		Vector<T>* ca_;
		uint32 index_;

	public:

		inline iterator(Vector<T>* ca, uint32 index) : ca_(ca), index_(index)
		{}
		inline iterator(const iterator& it) : ca_(it.ca_), index_(it.index_)
		{}
		inline iterator& operator=(const iterator& it)
		{
			ca_ = it.ca_;
			index_ = it.index_;
			return *this;
		}
		inline bool operator!=(iterator it) const
		{
			cgogn_assert(ca_ == it.ca_);
			return index_ != it.index_;
		}
		inline iterator& operator++()
		{
			index_ = ca_->container_->next_index(index_);
			return *this;
		}
		inline T& operator*() const
		{
			return ca_->operator[](index_);
		}
		inline uint32 index() { return index_; }
	};
	inline iterator begin() { return iterator(this, this->container_->first_index()); }
	inline iterator end() { return iterator(this, this->container_->last_index()); }
};

} // namespace cgogn

#endif // CGOGN_CORE_CONTAINER_VECTOR_H_

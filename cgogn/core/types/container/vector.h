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

#ifndef CGOGN_CORE_CONTAINER_VECTOR_H_
#define CGOGN_CORE_CONTAINER_VECTOR_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/numerics.h>

#include <cgogn/core/types/container/attribute_container.h>

#include <memory>
#include <string>
#include <vector>

namespace cgogn
{

//////////////////
// Vector class //
//////////////////

template <typename T>
class CGOGN_CORE_EXPORT Vector : public AttributeGenT
{
	using AttributeContainer = AttributeContainerT<Vector>;

private:
	std::vector<T> data_;

	inline void manage_index(uint32 index) override
	{
		while (index >= uint32(data_.size()))
			data_.push_back(T());
	}

public:
	Vector(AttributeContainer* container, const std::string& name) : AttributeGenT(container, name)
	{
		data_.reserve(512u);
	}

	~Vector() override
	{
	}

	inline T& operator[](uint32 index)
	{
		cgogn_message_assert(index < uint32(data_.size()), "index out of bounds");
		return data_[index];
	}

	inline const T& operator[](uint32 index) const
	{
		cgogn_message_assert(index < uint32(data_.size()), "index out of bounds");
		return data_[index];
	}

	inline void fill(const T& value)
	{
		std::fill(data_.begin(), data_.end(), value);
	}

	inline void swap(Vector<T>* vector)
	{
		if (vector->container_ == this->container_) // only swap from same container
			data_.swap(vector->data_);
	}

	inline void copy(Vector<T>* vector)
	{
		if (vector->container_ == this->container_) // only copy from same container
			data_ = vector->data_;
	}

	inline void clear() override
	{
		data_.clear();
		data_.reserve(512u);
		data_.shrink_to_fit();
	}

	inline std::shared_ptr<AttributeGenT> create_in(AttributeContainerGen& dst) const override
	{
		AttributeContainer* dst_container = dynamic_cast<AttributeContainer*>(&dst);
		if (dst_container)
		{
			auto attribute = dst_container->get_attribute<T>(name_);
			if (!attribute)
				attribute = dst_container->add_attribute<T>(name_);
			return attribute;
		}
		return nullptr;
	}

	inline void copy(const AttributeGenT& src) override
	{
		const Vector<T>* src_vector = dynamic_cast<const Vector<T>*>(&src);
		if (src_vector)
		{
			cgogn_message_assert(src_vector->data_.size() == data_.size(), "Copy from src with different capacity");
			data_ = src_vector->data_;
		}
	}

	inline const void* data_pointer() const
	{
		return &data_[0];
	}

	class const_iterator
	{
		const Vector<T>* vector_;
		uint32 index_;

	public:
		inline const_iterator(const Vector<T>* vector, uint32 index) : vector_(vector), index_(index)
		{
		}
		inline const_iterator(const const_iterator& it) : vector_(it.vector_), index_(it.index_)
		{
		}
		inline const_iterator& operator=(const const_iterator& it)
		{
			vector_ = it.vector_;
			index_ = it.index_;
			return *this;
		}
		inline bool operator!=(const_iterator it) const
		{
			cgogn_assert(vector_ == it.vector_);
			return index_ != it.index_;
		}
		inline const_iterator& operator++()
		{
			index_ = vector_->container_->next_index(index_);
			return *this;
		}
		inline const T& operator*() const
		{
			return vector_->operator[](index_);
		}
		inline uint32 index()
		{
			return index_;
		}
	};
	inline const_iterator begin() const
	{
		return const_iterator(this, this->container_->first_index());
	}
	inline const_iterator end() const
	{
		return const_iterator(this, this->container_->last_index());
	}

	class iterator
	{
		Vector<T>* vector_;
		uint32 index_;

	public:
		inline iterator(Vector<T>* vector, uint32 index) : vector_(vector), index_(index)
		{
		}
		inline iterator(const iterator& it) : vector_(it.vector_), index_(it.index_)
		{
		}
		inline iterator& operator=(const iterator& it)
		{
			vector_ = it.vector_;
			index_ = it.index_;
			return *this;
		}
		inline bool operator!=(iterator it) const
		{
			cgogn_assert(vector_ == it.vector_);
			return index_ != it.index_;
		}
		inline iterator& operator++()
		{
			index_ = vector_->container_->next_index(index_);
			return *this;
		}
		inline T& operator*() const
		{
			return vector_->operator[](index_);
		}
		inline uint32 index()
		{
			return index_;
		}
	};
	inline iterator begin()
	{
		return iterator(this, this->container_->first_index());
	}
	inline iterator end()
	{
		return iterator(this, this->container_->last_index());
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_CONTAINER_VECTOR_H_

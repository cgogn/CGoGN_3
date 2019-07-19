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

#ifndef CGOGN_CORE_CONTAINER_CHUNK_ARRAY_H_
#define CGOGN_CORE_CONTAINER_CHUNK_ARRAY_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/numerics.h>

#include <cgogn/core/types/container/attribute_container.h>

#include <vector>
#include <string>
#include <memory>

namespace cgogn
{

//////////////////////
// ChunkArray class //
//////////////////////

template <typename T>
class CGOGN_CORE_EXPORT ChunkArray : public AttributeGen
{
public:

	static const uint32 CHUNK_SIZE = 512;

private:

	std::vector<T*> chunks_;

	inline void manage_index(uint32 index) override
	{
		uint32 capacity = chunks_.size() * CHUNK_SIZE;
		while (index >= capacity)
		{
			chunks_.push_back(new T[CHUNK_SIZE]());
			capacity = chunks_.size() * CHUNK_SIZE;
		}
	}

public:

	ChunkArray(AttributeContainerGen* container, bool is_mark, const std::string& name) : AttributeGen(container, is_mark, name)
	{
		chunks_.reserve(512u);
	}

	~ChunkArray() override
	{
		for (auto chunk : chunks_)
			delete[] chunk;
	}

	inline T& operator[](uint32 index)
	{
		cgogn_message_assert(index / CHUNK_SIZE < chunks_.size(), "index out of bounds");
		return chunks_[index / CHUNK_SIZE][index % CHUNK_SIZE];
	}

	inline const T& operator[](uint32 index) const
	{
		cgogn_message_assert(index / CHUNK_SIZE < chunks_.size(), "index out of bounds");
		return chunks_[index / CHUNK_SIZE][index % CHUNK_SIZE];
	}

	inline void fill(const T& value)
	{
		for (auto chunk : chunks_)
			std::fill(chunk, chunk + CHUNK_SIZE, value);
	}

	inline void swap(ChunkArray<T>* ca)
	{
		if (ca->container_ == this->container_)
			chunks_.swap(ca->chunks_);
	}

	inline uint32 nb_chunks() const
	{
		return chunks_.size();
	}

	inline std::vector<const void*> chunk_pointers() const
	{
		std::vector<const void*> pointers;
		pointers.reserve(chunks_.size());
		for (auto chunk : chunks_)
			pointers.push_back(chunk);

		return pointers;
	}

	class const_iterator
	{
		const ChunkArray<T>* ca_;
		uint32 index_;

	public:

		inline const_iterator(const ChunkArray<T>* ca, uint32 index) : ca_(ca), index_(index)
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
		ChunkArray<T>* ca_;
		uint32 index_;

	public:

		inline iterator(ChunkArray<T>* ca, uint32 index) : ca_(ca), index_(index)
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

#endif // CGOGN_CORE_CONTAINER_CHUNK_ARRAY_H_

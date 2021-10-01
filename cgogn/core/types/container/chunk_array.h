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

#ifndef CGOGN_CORE_CONTAINER_CHUNK_ARRAY_H_
#define CGOGN_CORE_CONTAINER_CHUNK_ARRAY_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/numerics.h>

#include <cgogn/core/types/container/attribute_container.h>

#include <memory>
#include <string>
#include <vector>

namespace cgogn
{

//////////////////////
// ChunkArray class //
//////////////////////

template <typename T>
class CGOGN_CORE_EXPORT ChunkArray : public AttributeGenT
{
	using AttributeContainer = AttributeContainerT<ChunkArray>;

public:
	static const uint32 CHUNK_SIZE = 1024u;

private:
	std::vector<T*> chunks_;
	uint32 capacity_;

	inline void manage_index(uint32 index) override
	{
		while (index >= capacity_)
		{
			chunks_.push_back(new T[CHUNK_SIZE]());
			capacity_ = uint32(chunks_.size()) * CHUNK_SIZE;
		}
	}

public:
	ChunkArray(AttributeContainer* container, const std::string& name) : AttributeGenT(container, name)
	{
		chunks_.reserve(512u);
		capacity_ = 0u;
	}

	~ChunkArray() override
	{
		for (auto chunk : chunks_)
			delete[] chunk;
	}

	inline T& operator[](uint32 index)
	{
		cgogn_message_assert(index < capacity_, "index out of bounds");
		return chunks_[index / CHUNK_SIZE][index % CHUNK_SIZE];
	}

	inline const T& operator[](uint32 index) const
	{
		cgogn_message_assert(index < capacity_, "index out of bounds");
		return chunks_[index / CHUNK_SIZE][index % CHUNK_SIZE];
	}

	inline void fill(const T& value)
	{
		for (auto chunk : chunks_)
			std::fill(chunk, chunk + CHUNK_SIZE, value);
	}

	inline void swap(ChunkArray<T>* ca)
	{
		if (ca->container_ == this->container_) // only swap from same container
			chunks_.swap(ca->chunks_);
	}

	inline void copy(ChunkArray<T>* ca)
	{
		if (ca->container_ == this->container_) // only copy from same container
			for (uint32 i = 0; i < uint32(chunks_.size()); ++i)
				std::copy(ca->chunks_[i], ca->chunks_[i] + CHUNK_SIZE, chunks_[i]);
	}

	inline void clear() override
	{
		for (auto chunk : chunks_)
			delete[] chunk;
		chunks_.clear();
		capacity_ = 0;
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
		const ChunkArray<T>* src_ca = dynamic_cast<const ChunkArray<T>*>(&src);
		if (src_ca)
		{
			cgogn_message_assert(src_ca->capacity_ == capacity_, "Copy from src with different capacity");
			for (uint32 i = 0; i < uint32(src_ca->chunks_.size()); ++i)
				std::copy(src_ca->chunks_[i], src_ca->chunks_[i] + CHUNK_SIZE, chunks_[i]);
		}
	}

	inline uint32 nb_chunks() const
	{
		return uint32(chunks_.size());
	}

	inline std::vector<const void*> chunk_pointers() const
	{
		std::vector<const void*> pointers;
		pointers.reserve(uint32(chunks_.size()));
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
		{
		}
		inline const_iterator(const const_iterator& it) : ca_(it.ca_), index_(it.index_)
		{
		}
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
		ChunkArray<T>* ca_;
		uint32 index_;

	public:
		inline iterator(ChunkArray<T>* ca, uint32 index) : ca_(ca), index_(index)
		{
		}
		inline iterator(const iterator& it) : ca_(it.ca_), index_(it.index_)
		{
		}
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

#endif // CGOGN_CORE_CONTAINER_CHUNK_ARRAY_H_

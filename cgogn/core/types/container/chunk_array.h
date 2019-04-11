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
#include <cgogn/core/utils/assert.h>

#include <vector>
#include <string>
#include <memory>

namespace cgogn
{

template <uint32 CHUNK_SIZE>
class ChunkArrayContainerImpl;

template <uint32 CHUNK_SIZE>
class CGOGN_CORE_EXPORT ChunkArrayGenImpl
{
protected:

	std::string name_;
	ChunkArrayContainerImpl<CHUNK_SIZE>* container_;
	bool is_mark_;

	friend ChunkArrayContainerImpl<CHUNK_SIZE>;

	virtual void add_chunk() = 0;

public:

	ChunkArrayGenImpl(ChunkArrayContainerImpl<CHUNK_SIZE>* container, bool is_mark, const std::string& name) :
		name_(name),
		container_(container),
		is_mark_(is_mark)
	{}

	virtual ~ChunkArrayGenImpl()
	{
		if (is_mark_)
		{
			auto iter = std::find(container_->mark_chunk_arrays_.begin(), container_->mark_chunk_arrays_.end(), this);
			if (iter != container_->mark_chunk_arrays_.end())
			{
				*iter = container_->mark_chunk_arrays_.back();
				container_->mark_chunk_arrays_.pop_back();
			}
		}
		else
		{
			auto iter = std::find(container_->chunk_arrays_.begin(), container_->chunk_arrays_.end(), this);
			if (iter != container_->chunk_arrays_.end())
			{
				*iter = container_->chunk_arrays_.back();
				container_->chunk_arrays_.pop_back();
			}
		}
	}

	const std::string& name() const { return name_; }
	uint32 size() const { return container_->size(); }
	uint32 capacity() const { return container_->capacity(); }
};

template <typename T, uint32 CHUNK_SIZE>
class CGOGN_CORE_EXPORT ChunkArrayImpl : public ChunkArrayGenImpl<CHUNK_SIZE>
{
	std::vector<T*> chunks_;

	friend ChunkArrayContainerImpl<CHUNK_SIZE>;

	void add_chunk() override { chunks_.push_back(new T[CHUNK_SIZE]()); }

	void set_nb_chunks(uint32 nb_chunks)
	{
		if (nb_chunks >= chunks_.size())
		{
			for (std::size_t i = chunks_.size(); i < nb_chunks; ++i)
				add_chunk();
		}
		else
		{
			for (std::size_t i = static_cast<std::size_t>(nb_chunks); i < chunks_.size(); ++i)
				delete[] chunks_[i];
			chunks_.resize(nb_chunks);
		}
	}

public:

	class const_iterator
	{
		const ChunkArrayImpl<T, CHUNK_SIZE>* ca_;
		uint32 index_;

	public:

		inline const_iterator(const ChunkArrayImpl<T, CHUNK_SIZE>* ca, uint32 index) : ca_(ca), index_(index)
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
			++index_;
			return *this;
		}
		inline const T& operator*() const
		{
			return ca_->operator[](index_);
		}
		inline uint32 index() { return index_; }
	};
	inline const_iterator begin() const { return const_iterator(this, 0); }
	inline const_iterator end() const { return const_iterator(this, this->container_->size()); }

	class iterator
	{
		ChunkArrayImpl<T, CHUNK_SIZE>* ca_;
		uint32 index_;

	public:

		inline iterator(ChunkArrayImpl<T, CHUNK_SIZE>* ca, uint32 index) : ca_(ca), index_(index)
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
			++index_;
			return *this;
		}
		inline T& operator*() const
		{
			return ca_->operator[](index_);
		}
		inline uint32 index() { return index_; }
	};
	inline iterator begin() { return iterator(this, 0); }
	inline iterator end() { return iterator(this, this->container_->size()); }

	ChunkArrayImpl(ChunkArrayContainerImpl<CHUNK_SIZE>* container, bool is_mark, const std::string& name) : ChunkArrayGenImpl<CHUNK_SIZE>(container, is_mark, name)
	{
		chunks_.reserve(512u);
	}

	~ChunkArrayImpl() override
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

	inline void swap(std::shared_ptr<ChunkArrayImpl<T, CHUNK_SIZE>> ca)
	{
		if (ca->container_ == this->container_)
			chunks_.swap(ca->chunks_);
	}

	inline std::vector<const void*> chunk_pointers(uint32& chunk_byte_size) const
	{
		chunk_byte_size = CHUNK_SIZE * sizeof(T);

		std::vector<const void*> pointers;
		pointers.reserve(chunks_.size());
		for (auto c : chunks_)
			pointers.push_back(c);

		return pointers;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_CONTAINER_CHUNK_ARRAY_H_

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

#ifndef CGOGN_CORE_CONTAINER_CHUNK_ARRAY_CONTAINER_H_
#define CGOGN_CORE_CONTAINER_CHUNK_ARRAY_CONTAINER_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/numerics.h>

#include <vector>
#include <string>
#include <memory>

namespace cgogn
{

template <uint32 CHUNK_SIZE>
class CGOGN_CORE_EXPORT ChunkArrayContainer
{
public:

	/////////////////////////
	// ChunkArrayGen class //
	/////////////////////////

	class CGOGN_CORE_EXPORT ChunkArrayGen
	{
	protected:

		std::string name_;
		ChunkArrayContainer<CHUNK_SIZE>* container_;
		bool is_mark_;

		friend ChunkArrayContainer<CHUNK_SIZE>;

		virtual void add_chunk() = 0;

	public:

		ChunkArrayGen(ChunkArrayContainer<CHUNK_SIZE>* container, bool is_mark, const std::string& name) :
			name_(name),
			container_(container),
			is_mark_(is_mark)
		{}

		virtual ~ChunkArrayGen()
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

	//////////////////////
	// ChunkArray class //
	//////////////////////

	template <typename T>
	class CGOGN_CORE_EXPORT ChunkArray : public ChunkArrayGen
	{
		std::vector<T*> chunks_;

		friend ChunkArrayContainer<CHUNK_SIZE>;

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

		ChunkArray(ChunkArrayContainer<CHUNK_SIZE>* container, bool is_mark, const std::string& name) : ChunkArrayGen(container, is_mark, name)
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

		inline void swap(std::shared_ptr<ChunkArray<T>> ca)
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

	template <typename T>
	using ChunkArrayPtr = std::shared_ptr<ChunkArray<T>>;

	using ChunkArrayGenPtr = std::shared_ptr<ChunkArrayGen>;

private:

	std::vector<ChunkArrayGen*> chunk_arrays_;
	std::vector<ChunkArrayGenPtr> chunk_arrays_shared_ptr_;

	std::vector<ChunkArray<uint8>*> mark_chunk_arrays_;

	uint32 size_;
	uint32 nb_chunks_;

public:

	using const_iterator = typename std::vector<ChunkArrayGenPtr>::const_iterator;
	inline const_iterator begin() const { return chunk_arrays_shared_ptr_.begin(); }
	inline const_iterator end() const { return chunk_arrays_shared_ptr_.end(); }

	ChunkArrayContainer() : size_(0), nb_chunks_(0)
	{
		chunk_arrays_.reserve(16);
		chunk_arrays_shared_ptr_.reserve(16);
		mark_chunk_arrays_.reserve(16);
	}

	~ChunkArrayContainer()
	{}

	uint32 size() const { return size_; }
	uint32 capacity() const { return nb_chunks_ * CHUNK_SIZE; }

	template <typename T>
	ChunkArrayPtr<T> add_chunk_array(const std::string& name)
	{
		ChunkArrayPtr<T> casp = std::make_shared<ChunkArray<T>>(this, false, name);
		ChunkArray<T>* cap = casp.get();
		cap->set_nb_chunks(nb_chunks_);
		chunk_arrays_.push_back(cap);
		chunk_arrays_shared_ptr_.push_back(casp);
		return casp;
	}

	template <typename T>
	ChunkArrayPtr<T> get_chunk_array(const std::string& name) const
	{
		auto it = std::find_if(
			chunk_arrays_shared_ptr_.begin(),
			chunk_arrays_shared_ptr_.end(),
			[&] (const ChunkArrayGenPtr& att) { return att->name().compare(name) == 0; }
		);
		if (it != chunk_arrays_shared_ptr_.end())
			return std::dynamic_pointer_cast<ChunkArray<T>>(*it);
		return ChunkArrayPtr<T>();
	}

	void remove_chunk_array(ChunkArrayGenPtr attribute)
	{
		const std::string& name = attribute->name();
		auto it = std::find_if(
			chunk_arrays_shared_ptr_.begin(),
			chunk_arrays_shared_ptr_.end(),
			[&] (const ChunkArrayGenPtr& att) { return att->name().compare(name) == 0; }
		);
		if (it != chunk_arrays_shared_ptr_.end())
		{
			*it = chunk_arrays_shared_ptr_.back();
			chunk_arrays_shared_ptr_.pop_back();
		}
	}

	ChunkArray<uint8>* add_mark_chunk_array()
	{
		ChunkArray<uint8>* ca = new ChunkArray<uint8>(this, true, "mark");
		ca->set_nb_chunks(nb_chunks_);
		mark_chunk_arrays_.push_back(ca);
		return ca;
	}

	uint32 get_index()
	{
		uint32 index = size_++;
		if (index % CHUNK_SIZE == 0)
		{
			for (ChunkArrayGen* ag : chunk_arrays_)
				ag->add_chunk();
			for (ChunkArray<uint8>* a : mark_chunk_arrays_)
				a->add_chunk();
			++nb_chunks_;
		}
		return index;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_CONTAINER_CHUNK_ARRAY_CONTAINER_H_

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

#include <cgogn/core/types/container/chunk_array.h>
#include <cgogn/core/utils/numerics.h>

#include <vector>
#include <string>
#include <memory>

namespace cgogn
{

template <uint32 CHUNK_SIZE>
class CGOGN_CORE_EXPORT ChunkArrayContainerImpl
{
public:

	using ChunkArrayGen = ChunkArrayGenImpl<CHUNK_SIZE>;
	template <typename T>
	using ChunkArray = ChunkArrayImpl<T, CHUNK_SIZE>;

private:

	std::vector<ChunkArrayGen*> chunk_arrays_;
	std::vector<std::shared_ptr<ChunkArrayGen>> chunk_arrays_shared_ptr_;

	std::vector<ChunkArray<uint8>*> mark_chunk_arrays_;

	uint32 size_;
	uint32 nb_chunks_;

	friend ChunkArrayGen;

public:

	using const_iterator = typename std::vector<std::shared_ptr<ChunkArrayGen>>::const_iterator;
	inline const_iterator begin() const { return chunk_arrays_shared_ptr_.begin(); }
	inline const_iterator end() const { return chunk_arrays_shared_ptr_.end(); }

	ChunkArrayContainerImpl() : size_(0), nb_chunks_(0)
	{
		chunk_arrays_.reserve(16);
		chunk_arrays_shared_ptr_.reserve(16);
		mark_chunk_arrays_.reserve(16);
	}

	~ChunkArrayContainerImpl()
	{}

	uint32 size() const { return size_; }
	uint32 capacity() const { return nb_chunks_ * CHUNK_SIZE; }

	template <typename T>
	std::shared_ptr<ChunkArray<T>> add_chunk_array(const std::string& name)
	{
		std::shared_ptr<ChunkArray<T>> casp = std::make_shared<ChunkArray<T>>(this, false, name);
		ChunkArray<T>* cap = casp.get();
		cap->set_nb_chunks(nb_chunks_);
		chunk_arrays_.push_back(cap);
		chunk_arrays_shared_ptr_.push_back(casp);
		return casp;
	}

	template <typename T>
	std::shared_ptr<ChunkArray<T>> get_chunk_array(const std::string& name) const
	{
		auto it = std::find_if(
			chunk_arrays_shared_ptr_.begin(),
			chunk_arrays_shared_ptr_.end(),
			[&] (const std::shared_ptr<ChunkArrayGen>& att) { return att->name().compare(name) == 0; }
		);
		if (it != chunk_arrays_shared_ptr_.end())
			return std::dynamic_pointer_cast<ChunkArray<T>>(*it);
		return std::shared_ptr<ChunkArray<T>>();
	}

	void remove_chunk_array(std::shared_ptr<ChunkArrayGen> attribute)
	{
		const std::string& name = attribute->name();
		auto it = std::find_if(
			chunk_arrays_shared_ptr_.begin(),
			chunk_arrays_shared_ptr_.end(),
			[&] (const std::shared_ptr<ChunkArrayGen>& att) { return att->name().compare(name) == 0; }
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

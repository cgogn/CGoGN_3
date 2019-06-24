
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

#ifndef CGOGN_RENDERING_SHADERS_VBO_H_
#define CGOGN_RENDERING_SHADERS_VBO_H_

#include <GL/gl3w.h>
#include <string>
#include <iostream>
#include <vector>
#include <array>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/rendering_pureGL/vbo.h>

namespace cgogn
{

namespace rendering_pgl
{

/////////////////
// std::vector //
/////////////////

/**
  * @brief update vbo from a std::vector<VEC3>
  * @param vector
  * @param vbo vbo to update
  */
template <typename VEC3>
void update_vbo(const std::vector<VEC3>& vector, VBO* vbo)
{
	uint32 nb_elements = uint32(vector.size());
	vbo->allocate(nb_elements, 3);
	const uint32 vbo_bytes = nb_elements * 3 * uint32(sizeof(float32));

	// copy data
	vbo->bind();
	vbo->copy_data(0, vbo_bytes, vector.data());
	vbo->release();
}
template <typename VEC3, typename FUNC>
void update_vbo(const std::vector<VEC3>& vector, VBO* vbo, const FUNC& convert)
{
	using Vec3f = std::array<float32, 3>;
	uint32 nb_elements = uint32(vector.size());
	std::vector<Vec3f> converted;
	converted.reserve(nb_elements);
	for(const auto& v: vector)
		converted.push_back(convert(v));

	vbo->allocate(nb_elements, 3);
	const uint32 vbo_bytes = nb_elements * 3 * uint32(sizeof(float32));

	// copy data
	vbo->bind();
	vbo->copy_data(0, vbo_bytes, converted.data());
	vbo->release();
}

////////////
// Vector //
////////////

/**
 * @brief update vbo from a Vector<VEC3>
 * @param attribute
 * @param vbo vbo to update
 */
template <typename VEC3>
void update_vbo(const Vector<VEC3>* attribute, VBO* vbo)
{
	vbo->set_name(attribute->name());

	uint32 nb_elements = attribute->maximum_index();

	vbo->allocate(nb_elements, 3);

	// copy data
	vbo->bind();
	vbo->copy_data(0, nb_elements * sizeof (VEC3), attribute->data_pointer());
	vbo->release();
}

////////////////
// ChunkArray //
////////////////

/**
 * @brief update vbo from a Vector<VEC3>
 * @param attribute
 * @param vbo vbo to update
 */
template <typename VEC3>
void update_vbo(const ChunkArray<VEC3>* attribute, VBO* vbo)
{
	vbo->set_name(attribute->name());

	uint32 nb_elements = attribute->maximum_index();

	uint32 chunk_size = ChunkArray<VEC3>::CHUNK_SIZE;
	vbo->allocate(nb_elements + (chunk_size - nb_elements % chunk_size), 3);

	uint32 chunk_byte_size;
	std::vector<const void*> chunk_pointers = attribute->chunk_pointers(chunk_byte_size);

	// copy data
	vbo->bind();
	for (uint32 i = 0, size = uint32(chunk_pointers.size()); i < size; ++i)
		vbo->copy_data(i * chunk_byte_size, chunk_byte_size, chunk_pointers[i]);
	vbo->release();
}


template <typename VEC3, typename FUNC>
void update_vbo(const ChunkArray<VEC3>* attribute, VBO* vbo, const FUNC& convert)
{
	vbo->set_name(attribute->name());

	uint32 nb_elements = attribute->maximum_index();

	uint32 chunk_size = ChunkArray<VEC3>::CHUNK_SIZE;
	uint32 chunk_byte_size;
	std::vector<const void*> chunk_pointers = attribute->chunk_pointers(chunk_byte_size);


	// copy (after conversion)
	using OutputConvert = func_return_type<FUNC>;
	vbo->allocate(nb_elements + (chunk_size - nb_elements % chunk_size), sizeof(OutputConvert)/4);

	OutputConvert* dst = reinterpret_cast<OutputConvert*>(vbo->lock_pointer());
	for (uint32 i = 0; i < chunk_pointers.size(); ++i)
	{
		const VEC3* typed_chunk = static_cast<const VEC3*>(chunk_pointers[i]);
		for (uint32 j = 0; j < chunk_size; ++j)
			*dst++ = convert(typed_chunk[j]);
	}
	vbo->release_pointer();
}



} // namespace rendering_pgl

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_VBO_H_


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

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/rendering/vbo.h>

#include <GL/gl3w.h>

#include <string>
#include <iostream>
#include <vector>
#include <array>

namespace cgogn
{

namespace rendering
{

/////////////////
// std::vector //
/////////////////

/**
  * @brief update vbo from a std::vector<VEC>
  * @param vector
  * @param vbo vbo to update
  */
template <typename VEC,
		  typename std::enable_if<std::is_same<typename geometry::vector_traits<VEC>::Scalar, float32>::value>::type* = nullptr>
void update_vbo(const std::vector<VEC>& vector, VBO* vbo)
{
	static const std::size_t element_size = geometry::vector_traits<VEC>::SIZE;
	uint32 nb_elements = uint32(vector.size());
	vbo->allocate(nb_elements, element_size);

	// copy data
	vbo->bind();
	vbo->copy_data(0, nb_elements * element_size * uint32(sizeof(float32)), vector.data());
	vbo->release();
}

/**
 * @brief update vbo from an on-the-fly converted std::vector<VEC>
 * @param vector
 * @param vbo vbo to update
 * @param convert the conversion function
 */
template <typename VEC, typename FUNC>
void update_vbo(const std::vector<VEC>& vector, VBO* vbo, const FUNC& convert)
{
	static_assert(is_func_parameter_same<FUNC, const VEC&>::value, "Wrong conversion function parameter type");
	
	using OutputType = func_return_type<FUNC>;
	static const std::size_t output_type_size = geometry::vector_traits<OutputType>::SIZE;

	uint32 nb_elements = uint32(vector.size());
	vbo->allocate(nb_elements, output_type_size);

	OutputType* dst = reinterpret_cast<OutputType*>(vbo->lock_pointer());
	for (const VEC& v : vector)
		*dst++ = convert(v);
	vbo->release_pointer();
}

template <typename VEC,
		  typename std::enable_if<std::is_same<typename geometry::vector_traits<VEC>::Scalar, float64>::value>::type* = nullptr>
void update_vbo(const std::vector<VEC>* vector, VBO* vbo)
{
	static const std::size_t element_size = geometry::vector_traits<VEC>::SIZE;
	switch (element_size)
	{
		case 1:
			update_vbo(vector, vbo, [] (const VEC& n) -> float32
			{
				return float32(n[0]);
			});
			break;
		case 2:
			update_vbo(vector, vbo, [] (const VEC& n) -> geometry::Vec2f
			{
				return { float32(n[0]), float32(n[1]) };
			});
			break;
		case 3:
			update_vbo(vector, vbo, [] (const VEC& n) -> geometry::Vec3f
			{
				return { float32(n[0]), float32(n[1]), float32(n[2]) };
			});
			break;
		case 4:
			update_vbo(vector, vbo, [] (const VEC& n) -> geometry::Vec4f
			{
				return { float32(n[0]), float32(n[1]), float32(n[2]), float32(n[3]) };
			});
			break;
	}
}

////////////
// Vector //
////////////

/**
 * @brief update vbo from a Vector<VEC>
 * @param attribute
 * @param vbo vbo to update
 */
template <typename VEC,
		  typename std::enable_if<std::is_same<typename geometry::vector_traits<VEC>::Scalar, float32>::value>::type* = nullptr>
void update_vbo(const Vector<VEC>* attribute, VBO* vbo)
{
	vbo->set_name(attribute->name());

	static const std::size_t element_size = geometry::vector_traits<VEC>::SIZE;
	uint32 nb_elements = attribute->maximum_index();
	vbo->allocate(nb_elements, element_size);

	// copy data
	vbo->bind();
	vbo->copy_data(0, nb_elements * element_size * uint32(sizeof(float32)), attribute->data_pointer());
	vbo->release();
}

/**
 * @brief update vbo from an on-the-fly converted Vector<VEC>
 * @param attribute
 * @param vbo vbo to update
 * @param convert the conversion function
 */
template <typename VEC, typename FUNC>
void update_vbo(const Vector<VEC>* attribute, VBO* vbo, const FUNC& convert)
{
	static_assert(is_func_parameter_same<FUNC, const VEC&>::value, "Wrong conversion function parameter type");
	
	vbo->set_name(attribute->name());

	using OutputType = func_return_type<FUNC>;
	static const std::size_t output_type_size = geometry::vector_traits<OutputType>::SIZE;

	uint32 nb_elements = attribute->maximum_index();
	vbo->allocate(nb_elements, output_type_size);

	OutputType* dst = reinterpret_cast<OutputType*>(vbo->lock_pointer());
	for (const VEC& v : *attribute)
		*dst++ = convert(v);
	vbo->release_pointer();
}

template <typename VEC,
		  typename std::enable_if<std::is_same<typename geometry::vector_traits<VEC>::Scalar, float64>::value>::type* = nullptr>
void update_vbo(const Vector<VEC>* attribute, VBO* vbo)
{
	static const std::size_t element_size = geometry::vector_traits<VEC>::SIZE;
	switch (element_size)
	{
		case 1:
			update_vbo(attribute, vbo, [] (const VEC& n) -> float32
			{
				return float32(n[0]);
			});
			break;
		case 2:
			update_vbo(attribute, vbo, [] (const VEC& n) -> geometry::Vec2f
			{
				return { float32(n[0]), float32(n[1]) };
			});
			break;
		case 3:
			update_vbo(attribute, vbo, [] (const VEC& n) -> geometry::Vec3f
			{
				return { float32(n[0]), float32(n[1]), float32(n[2]) };
			});
			break;
		case 4:
			update_vbo(attribute, vbo, [] (const VEC& n) -> geometry::Vec4f
			{
				return { float32(n[0]), float32(n[1]), float32(n[2]), float32(n[3]) };
			});
			break;
	}
}

////////////////
// ChunkArray //
////////////////

template <typename VEC,
		  typename std::enable_if<std::is_same<typename geometry::vector_traits<VEC>::Scalar, float32>::value>::type* = nullptr>
void update_vbo(const ChunkArray<VEC>* attribute, VBO* vbo)
{
	vbo->set_name(attribute->name());

	static const std::size_t element_size = geometry::vector_traits<VEC>::SIZE;
	static const uint32 chunk_size = ChunkArray<VEC>::CHUNK_SIZE;
	uint32 nb_chunks = attribute->nb_chunks();
	vbo->allocate(nb_chunks * chunk_size, element_size);

	std::vector<const void*> chunk_pointers = attribute->chunk_pointers();

	// copy data
	uint32 vbo_chunk_byte_size = chunk_size * element_size * uint32(sizeof(float32));
	vbo->bind();
	for (uint32 i = 0, size = uint32(chunk_pointers.size()); i < size; ++i)
		vbo->copy_data(i * vbo_chunk_byte_size, vbo_chunk_byte_size, chunk_pointers[i]);
	vbo->release();
}

/**
 * @brief update vbo from an on-the-fly converted ChunkArray<VEC>
 * @param attribute
 * @param vbo vbo to update
 * @param convert the conversion function
 */
template <typename VEC, typename FUNC>
void update_vbo(const ChunkArray<VEC>* attribute, VBO* vbo, const FUNC& convert)
{
	static_assert(is_func_parameter_same<FUNC, const VEC&>::value, "Wrong conversion function parameter type");
	
	vbo->set_name(attribute->name());

	using OutputType = func_return_type<FUNC>;
	static const std::size_t output_type_size = geometry::vector_traits<OutputType>::SIZE;
	static const uint32 chunk_size = ChunkArray<VEC>::CHUNK_SIZE;
	
	uint32 nb_elements = attribute->maximum_index();
	vbo->allocate(nb_elements, output_type_size);

	std::vector<const void*> chunk_pointers = attribute->chunk_pointers();

	OutputType* dst = reinterpret_cast<OutputType*>(vbo->lock_pointer());
	for (uint32 i = 0, size = uint32(chunk_pointers.size()); i < size; ++i)
	{
		const VEC* chunk = static_cast<const VEC*>(chunk_pointers[i]);
		for (uint32 j = 0; j < chunk_size && i * chunk_size + j < nb_elements; ++j)
			*dst++ = convert(chunk[j]);
	}
	vbo->release_pointer();
}

template <typename VEC,
		  typename std::enable_if<std::is_same<typename geometry::vector_traits<VEC>::Scalar, float64>::value>::type* = nullptr>
void update_vbo(const ChunkArray<VEC>* attribute, VBO* vbo)
{
	static const std::size_t element_size = geometry::vector_traits<VEC>::SIZE;
	switch (element_size)
	{
		case 1:
			update_vbo(attribute, vbo, [] (const VEC& n) -> float32
			{
				return float32(n[0]);
			});
			break;
		case 2:
			update_vbo(attribute, vbo, [] (const VEC& n) -> geometry::Vec2f
			{
				return { float32(n[0]), float32(n[1]) };
			});
			break;
		case 3:
			update_vbo(attribute, vbo, [] (const VEC& n) -> geometry::Vec3f
			{
				return { float32(n[0]), float32(n[1]), float32(n[2]) };
			});
			break;
		case 4:
			update_vbo(attribute, vbo, [] (const VEC& n) -> geometry::Vec4f
			{
				return { float32(n[0]), float32(n[1]), float32(n[2]), float32(n[3]) };
			});
			break;
	}
}

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_VBO_H_


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

#include <QOpenGLBuffer>

#include <iostream>
#include <vector>
#include <array>

namespace cgogn
{

namespace rendering
{

class VBO
{
protected:

	std::size_t nb_vectors_;
	uint32 vector_dimension_;
	QOpenGLBuffer buffer_;
	std::string name_;

public:

	inline VBO(uint32 vec_dim = 3u) :
		nb_vectors_(),
		vector_dimension_(vec_dim)
	{
		const bool buffer_created = buffer_.create();
		if (!buffer_created)
		{
			std::cerr << "VBO::VBO(uint32): The call to QOpenGLBuffer::create() failed. Maybe there is no QOpenGLContext." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		buffer_.bind();
		buffer_.setUsagePattern(QOpenGLBuffer::StreamDraw);
	}

	inline ~VBO()
	{
		buffer_.destroy();
	}

	inline void set_name(const std::string& name) { name_ = name; }

	inline const std::string& name() const { return name_; }

	inline void bind() { buffer_.bind(); }

	inline void release() { buffer_.release(); }

	/**
	 * @brief allocate VBO memory
	 * @param nb_vectors number of vectors
	 * @param vector_dimension_ number of component of each vector
	 */
	inline void allocate(std::size_t nb_vectors, uint32 vector_dimension)
	{
		buffer_.bind();
		std::size_t total = nb_vectors * vector_dimension;
		if (total != nb_vectors_ * vector_dimension_) // only allocate when > ?
			buffer_.allocate(uint32(total) * uint32(sizeof(float32)));
		nb_vectors_ = nb_vectors;
		if (vector_dimension != vector_dimension_)
		{
			vector_dimension_ = vector_dimension;
			std::cerr << "VBO::allocate: Changing the VBO vector_dimension." << std::endl;
		}
		buffer_.release();
	}

	/**
	 * @brief get and lock pointer on buffer memory
	 * @return  the pointer
	 */
	inline float32* lock_pointer()
	{
		buffer_.bind();
		return reinterpret_cast<float32*>(buffer_.map(QOpenGLBuffer::ReadWrite));
	}

	/**
	 * @brief release_pointer
	 */
	inline void release_pointer() { buffer_.unmap(); }

	/**
	 * @brief copy data
	 * @param offset offset in bytes in the bufffer
	 * @param nb number of bytes to copy
	 * @param src source pointer
	 */
	inline void copy_data(uint32 offset, std::size_t nb, const void* src)
	{
		buffer_.write(int(offset), src, int(nb));
	}

	/**
	 * @brief dimension of vectors stored in buffer
	 */
	inline uint32 vector_dimension() const { return vector_dimension_; }

	uint32 size() const { return uint32(nb_vectors_); }

	GLuint id() const { return buffer_.bufferId(); }
};

////////////
// vector //
////////////

/**
  * @brief update vbo from a std::vector<VEC>
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

////////////////
// ChunkArray //
////////////////

/**
 * @brief update vbo from one Attribute
 * @param attribute Attribute (must contain Vec3<float>)
 * @param vbo vbo to update
 */
template <typename VEC3>
void update_vbo(CMapBase::AttributePtr<VEC3> attribute, VBO* vbo)
{
	vbo->set_name(attribute->name());

	uint32 nb_elements = attribute->maximum_index();

	// TODO: remove chunk array hack

	vbo->allocate(nb_elements + (512 - nb_elements % 512), 3);

	uint32 chunk_byte_size;
	std::vector<const void*> chunk_pointers = attribute->chunk_pointers(chunk_byte_size);

	// copy data
	vbo->bind();
	for (uint32 i = 0, size = uint32(chunk_pointers.size()); i < size; ++i)
		vbo->copy_data(i * chunk_byte_size, chunk_byte_size, chunk_pointers[i]);
	vbo->release();
}

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_VBO_H_

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
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  021binding_point-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/

#include <cgogn/rendering/mesh_render.h>

namespace cgogn
{

namespace rendering
{

MeshRender::MeshRender()
{
	for (uint32 i = 0u; i < SIZE_BUFFER; ++i)
	{
		indices_buffers_[i] = std::make_unique<EBO>();
		indices_buffers_uptodate_[i] = false;
	}
}

MeshRender::~MeshRender()
{
}

void MeshRender::draw(DrawingType prim, GLint binding_point)
{
	uint32 prim_buffer = (prim < SIZE_BUFFER) ? prim : prim - (SIZE_BUFFER + 1);
	int32 nb_indices = int32(indices_buffers_[prim_buffer]->size());
	if (nb_indices == 0)
		return;

	switch (prim)
	{
	case POINTS:
		indices_buffers_[prim]->bind();
		glDrawElements(GL_POINTS, nb_indices, GL_UNSIGNED_INT, nullptr);
		indices_buffers_[prim]->release();
		break;
	case LINES:
		indices_buffers_[prim]->bind();
		glDrawElements(GL_LINES, nb_indices, GL_UNSIGNED_INT, nullptr);
		indices_buffers_[prim]->release();
		break;
	case TRIANGLES:
		indices_buffers_[prim]->bind();
		glDrawElements(GL_TRIANGLES, nb_indices, GL_UNSIGNED_INT, nullptr);
		indices_buffers_[prim]->release();
		break;
	case VOLUMES_EDGES:
		indices_buffers_[prim]->bind_tb(binding_point);
		glDrawArraysInstanced(GL_LINES, 0, 2, nb_indices / 3);
		indices_buffers_[prim]->release_tb();
		break;
	case VOLUMES_FACES:
		indices_buffers_[prim]->bind_tb(binding_point);
		glDrawArraysInstanced(GL_TRIANGLES, 0, 3, nb_indices / 4);
		indices_buffers_[prim]->release_tb();
		break;
	case VOLUMES_VERTICES:
		indices_buffers_[prim]->bind_tb(binding_point);
		glDrawArrays(GL_POINTS, 0, nb_indices / 2);
		indices_buffers_[prim]->release_tb();
		break;
	case LINES_TB:
		indices_buffers_[LINES]->bind_tb(binding_point);
		indices_buffers_[INDEX_EDGES]->bind_tb(binding_point + 1);
		glDrawArraysInstanced(GL_LINES, 0, 2, nb_indices / 2);
		indices_buffers_[INDEX_EDGES]->release_tb();
		indices_buffers_[LINES]->release_tb();
		break;
	case TRIANGLES_TB:
		indices_buffers_[TRIANGLES]->bind_tb(binding_point);
		indices_buffers_[INDEX_FACES]->bind_tb(binding_point + 1);
		glDrawArraysInstanced(GL_TRIANGLES, 0, 3, nb_indices / 3);
		indices_buffers_[INDEX_FACES]->release_tb();
		indices_buffers_[TRIANGLES]->release_tb();
		break;
	default:
		break;
	}
}

} // namespace rendering

} // namespace cgogn

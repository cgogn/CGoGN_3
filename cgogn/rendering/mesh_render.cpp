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
		nb_indices_[i] = 0;
	}
}

MeshRender::~MeshRender()
{
}

void MeshRender::draw(DrawingBufferType prim, GLint binding_point)
{
	if (nb_indices_[prim] == 0)
		return;

	switch (prim)
	{
	case BUFFER_POINTS:
		indices_buffers_[prim]->bind();
		glDrawElements(GL_POINTS, nb_indices_[prim], GL_UNSIGNED_INT, nullptr);
		indices_buffers_[prim]->release();
		break;
	case BUFFER_LINES:
		indices_buffers_[prim]->bind();
		glDrawElements(GL_LINES, nb_indices_[prim], GL_UNSIGNED_INT, nullptr);
		indices_buffers_[prim]->release();
		break;
	case BUFFER_LINES_TB:
		indices_buffers_[BUFFER_LINES]->bind_tb(binding_point);
		glDrawArraysInstanced(GL_LINES, 0, 2, nb_indices_[BUFFER_LINES]/2);
		indices_buffers_[BUFFER_LINES]->release_tb();
		break;
	case BUFFER_TRIANGLES:
		indices_buffers_[prim]->bind();
		glDrawElements(GL_TRIANGLES, nb_indices_[prim], GL_UNSIGNED_INT, nullptr);
		indices_buffers_[prim]->release();
		break;
	case BUFFER_TRIANGLES_TB:
		indices_buffers_[BUFFER_TRIANGLES]->bind_tb(binding_point);
		glDrawArraysInstanced(GL_TRIANGLES, 0, 3, nb_indices_[BUFFER_TRIANGLES]/3);
		indices_buffers_[BUFFER_TRIANGLES]->release_tb();
		break;
	case BUFFER_VOLUMES_EDGES:
		indices_buffers_[prim]->bind_tb(binding_point);
		glDrawArraysInstanced(GL_LINES, 0, 2, nb_indices_[prim]/3);
		indices_buffers_[prim]->release_tb();
		break;
	case BUFFER_VOLUMES_FACES:
		indices_buffers_[prim]->bind_tb(binding_point);
		glDrawArraysInstanced(GL_TRIANGLES, 0, 3, nb_indices_[prim]/4);
		indices_buffers_[prim]->release_tb();
		break;
	case BUFFER_VOLUMES_VERTICES:
		indices_buffers_[prim]->bind_tb(binding_point);
		glDrawArrays(GL_POINTS, 0, nb_indices_[prim]/2);
		indices_buffers_[prim]->release_tb();
		break;
	default:
		break;
	}
	GL_ASSERT("")
}


} // namespace rendering

} // namespace cgogn

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

void MeshRender::draw(DrawingType prim)
{
	if (nb_indices_[prim] == 0)
		return;

	indices_buffers_[prim]->bind();
	switch (prim)
	{
	case POINTS:
		indices_buffers_[POINTS]->bind();
		glDrawElements(GL_POINTS, nb_indices_[POINTS], GL_UNSIGNED_INT, 0);
		indices_buffers_[POINTS]->release();
		break;
	case LINES:
		indices_buffers_[LINES]->bind();
		glDrawElements(GL_LINES, nb_indices_[LINES], GL_UNSIGNED_INT, 0);
		indices_buffers_[LINES]->release();
		break;
	case TRIANGLES:
		indices_buffers_[TRIANGLES]->bind();
		glDrawElements(GL_TRIANGLES, nb_indices_[TRIANGLES], GL_UNSIGNED_INT, 0);
		indices_buffers_[TRIANGLES]->release();
		break;
	case BOUNDARY:
		break;
	case VOLUMES_EDGES:
		indices_buffers_[VOLUMES_EDGES]->bind_tb(10);
		glDrawArraysInstanced(GL_LINES, 0, 2, nb_indices_[VOLUMES_EDGES]/3);
		indices_buffers_[VOLUMES_EDGES]->release();
		break;
	case VOLUMES_FACES:
		indices_buffers_[VOLUMES_FACES]->bind_tb(10);
		glDrawArraysInstanced(GL_TRIANGLES, 0, 3, nb_indices_[VOLUMES_FACES]/4);
		indices_buffers_[VOLUMES_FACES]->release();
		break;
	case VOLUMES_VERTICES:
		indices_buffers_[VOLUMES_VERTICES]->bind_tb(10);
		glDrawArrays(GL_POINTS, 0, nb_indices_[VOLUMES_VERTICES]/2);
		indices_buffers_[VOLUMES_VERTICES]->release();
		break;
	default:
		break;
	}
	indices_buffers_[prim]->release();
}

} // namespace rendering

} // namespace cgogn

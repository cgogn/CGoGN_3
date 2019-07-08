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

#include <cgogn/rendering_pureGL/mesh_render.h>

namespace cgogn
{

namespace rendering_pgl
{

MeshRender::MeshRender()
{
	for (uint32 i = 0u; i < SIZE_BUFFER; ++i)
	{
		indices_buffers_[i] = make_unique<EBO>();
		indices_buffers_uptodate_[i] = false;
		nb_indices_[i] = 0;
	}
}

MeshRender::~MeshRender()
{}

void MeshRender::draw(DrawingType prim)
{
	if (nb_indices_[prim] == 0)
		return;

	indices_buffers_[prim]->bind();
	switch (prim)
	{
		case POINTS:
			glDrawElements(GL_POINTS, nb_indices_[POINTS], GL_UNSIGNED_INT, 0);
			break;
		case LINES:
			glDrawElements(GL_LINES, nb_indices_[LINES], GL_UNSIGNED_INT, 0);
			break;
		case TRIANGLES:
			glDrawElements(GL_TRIANGLES, nb_indices_[TRIANGLES], GL_UNSIGNED_INT, 0);
			break;
		case BOUNDARY:
			default:
			break;
	}
	indices_buffers_[prim]->release();
}

} // namespace rendering_pgl

} // namespace cgogn

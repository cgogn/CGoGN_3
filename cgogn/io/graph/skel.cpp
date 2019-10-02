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

#include <cgogn/io/graph/cg.h>
#include <cgogn/io/utils.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>

#include <cgogn/geometry/types/vector_traits.h>

#include <vector>
#include <fstream>


#include <cgogn/core/types/cmap/cmap_info.h>

namespace cgogn
{

namespace io
{

bool import_SKEL(Graph& g, const std::string& filename)
{
	return false;
}

} // namespace io

} // namespace cgogn

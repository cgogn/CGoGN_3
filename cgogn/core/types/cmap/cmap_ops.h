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

#ifndef CGOGN_CORE_CMAP_CMAP_OPS_H_
#define CGOGN_CORE_CMAP_CMAP_OPS_H_

namespace cgogn
{

template <typename CELL, typename CMAP,
		  typename = typename std::enable_if<std::is_base_of<CMapBase, CMAP>::value>::type>
void
set_embedding(CMAP& m, CELL c, uint32 emb)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<CMAP>::Cells>::value, "CELL not supported in this MESH");
	m.foreach_dart_of_orbit(c, [&] (Dart d) -> bool
	{
		m.template set_embedding<CELL>(d, emb);
		return true;
	});
}

template <typename CELL, typename CMAP,
		  typename = typename std::enable_if<std::is_base_of<CMapBase, CMAP>::value>::type>
void
create_embedding(CMAP& m, CELL c)
{
	static_assert(is_in_tuple<CELL, typename mesh_traits<CMAP>::Cells>::value, "CELL not supported in this MESH");
	uint32 emb = m.attribute_containers_[CELL::ORBIT].new_index();
	set_embedding(m, c, emb);
}

}

#endif // CGOGN_CORE_CMAP_CMAP_OPS_H_

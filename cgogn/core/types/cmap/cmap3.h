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

#ifndef CGOGN_CORE_TYPES_CMAP_CMAP3_H_
#define CGOGN_CORE_TYPES_CMAP_CMAP3_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cmap2.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMap3 : public CMap2
{
	using Self = CMap3;
	using Inherit = CMap2;
	
	std::shared_ptr<Attribute<Dart>> phi3_;

	using Vertex = Cell<PHI21_PHI31>;
	using Vertex2 = Cell<PHI21>;
	using Edge = Cell<PHI2_PHI3>;
	using Edge2 = Inherit::Edge;
	using Face = Cell<PHI1_PHI3>;
	using Face2 = Cell<PHI1>;
	using Volume = Cell<PHI1_PHI2>;
	using CC = Cell<PHI1_PHI2_PHI3>;

	using Cells = std::tuple<Vertex, Vertex2, Edge, Edge2, Face, Face2, Volume>;

	CMap3() : CMap2()
	{
		phi3_ = base_map_->add_or_get_relation("phi3");
	}
	CMap3(std::shared_ptr<CMapBase> m) : CMap2(m)
	{
		phi3_ = base_map_->add_or_get_relation("phi3");
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CMAP3_H_

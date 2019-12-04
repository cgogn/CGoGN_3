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

#ifndef CGOGN_CORE_TYPES_CMAP_CMAP2_H_
#define CGOGN_CORE_TYPES_CMAP_CMAP2_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cmap1.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMap2 : public CMap1
{
	std::shared_ptr<Attribute<Dart>> phi2_;

	using Vertex = Cell<PHI21>;
	using HalfEdge = Cell<DART>;
	using Edge = Cell<PHI2>;
	using Face = Cell<PHI1>;
	using Volume = Cell<PHI1_PHI2>;
	using CC = Volume;

	using Cells = std::tuple<Vertex, HalfEdge, Edge, Face, Volume>;
	
	using AttributeContainer = Inherit::AttributeContainer;

	template <typename T>
	using Attribute = Inherit::Attribute<T>;
	using AttributeGen = Inherit::AttributeGen;
	using MarkAttribute = Inherit::MarkAttribute;

	CMap2() : CMap1()
	{
		phi2_ = base_map_->add_or_get_relation("phi2");
	}
	CMap2(std::shared_ptr<CMapBase> m) : CMap1(m)
	{
		phi2_ = base_map_->add_or_get_relation("phi2");
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CMAP2_H_

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

#ifndef CGOGN_CORE_TYPES_CMAP_CMAP0_H_
#define CGOGN_CORE_TYPES_CMAP_CMAP0_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/cmap_base.h>
#include <memory>

namespace cgogn
{

struct CGOGN_CORE_EXPORT CMap0
{
	using Self = CMap0;
	
	using Vertex = Cell<DART>;
	using CC = Vertex;

	using Cells = std::tuple<Vertex>;
	
	using AttributeContainer = CMapBase::AttributeContainer;

	template <typename T>
	using Attribute = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	using MarkAttribute = CMapBase::MarkAttribute;
	
	static const bool is_mesh_view = true;
	
   std::shared_ptr<CMapBase> base_map_;

    CMap0()
    {
        base_map_ = std::make_shared<CMapBase>();
    }

    CMap0(std::shared_ptr<CMapBase> m) : base_map_(m)
    {}
	
    inline const std::shared_ptr<CMapBase> mesh() const {return base_map_;}
    inline std::shared_ptr<CMapBase> mesh(){return base_map_;}
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CMAP0_H_

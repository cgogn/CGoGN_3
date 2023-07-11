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

#ifndef CGOGN_CORE_TYPES_GMAP_GMAP0_H_
#define CGOGN_CORE_TYPES_GMAP_GMAP0_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/maps/gmap/gmap_base.h>
//#include <cgogn/core/types/maps/cell.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT GMap0 : public GMapBase
{
	static const uint8 dimension = 0;

	using Vertex = Cell < Orbit::DART > ;
	using Edge = Cell < Orbit::BETA0 > ;

	using Cells = std::tuple<Vertex,Edge>;

	std::shared_ptr<Attribute<Dart>> beta0_;

	GMap0()
	{
		beta0_ = add_relation("beta0");
	}

};

template <>
struct mesh_traits<GMap0>
{
	using MeshType = GMap0;
	static constexpr const char* name = "GMap0";
	static constexpr const uint8 dimension = 0;

	using Vertex = typename GMap0::Vertex;
	using Edge = typename GMap0::Edge;

	using Cells = std::tuple<Vertex,Edge>;
	static constexpr const char* cell_names[] = {"Vertex","Edge"};

	template <typename T>
	using Attribute = MapBase::Attribute<T>;
	using AttributeGen = MapBase::AttributeGen;
	using MarkAttribute = MapBase::MarkAttribute;
};

inline Dart beta0(const GMap0& m, Dart d)
{
	return (*(m.beta0_))[d.index];
}

inline void beta0_sew(GMap0& m, Dart d, Dart e)
{
	cgogn_assert(beta0(m, d) == d);
	cgogn_assert(beta0(m, e) == e);
	(*(m.beta0_))[d.index] = e;
	(*(m.beta0_))[e.index] = d;
}

inline void beta0_unsew(GMap0& m, Dart d)
{
	Dart e = beta0(m, d);
	(*(m.beta0_))[d.index] = d;
	(*(m.beta0_))[e.index] = e;
}


GMap0::Edge add_edge(GMap0& m, bool set_indices = true);

void remove_edge(GMap0& m, GMap0::Edge e, bool set_indice = true);

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_GMAP_GMAP0_H_

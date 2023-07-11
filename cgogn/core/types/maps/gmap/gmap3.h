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

#ifndef CGOGN_CORE_TYPES_GMAP_GMAP3_H_
#define CGOGN_CORE_TYPES_GMAP_GMAP3_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/maps/gmap/gmap2.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT GMap3 : public GMap2
{
	static const uint8 dimension = 3;

	using Parent = GMap2;

	using Vertex = Cell<Orbit::BETA1_BETA2_BETA3>;
	using Vertex2 = Cell<Orbit::BETA1_BETA2>;
	using HalfEdge = Cell<Orbit::BETA0>;
	using Edge = Cell<Orbit::BETA0_BETA2_BETA3>;
	using Edge2 = Cell<Orbit::BETA0_BETA2>;
	using Face = Cell<Orbit::BETA0_BETA1_BETA3>;
	using Face2 = Cell<Orbit::BETA0_BETA1>;
	using Volume = Cell<Orbit::BETA0_BETA1_BETA2>;
	using CC = Cell<Orbit::BETA0_BETA1_BETA2_BETA3>;
	using Cells = std::tuple<Vertex, Vertex2, HalfEdge, Edge, Edge2, Face, Face2, Volume>;

	std::shared_ptr<Attribute<Dart>> beta3_;


	inline GMap3() : GMap2()
	{
		beta3_ = add_relation("beta3");
	}


};

template <>
struct mesh_traits<GMap3>
{
	using MeshType = GMap3;
	static constexpr const char* name = "GMap3";
	static constexpr const uint8 dimension = 3;

	using Parent = GMap3::Parent;

	using Vertex = GMap3::Vertex;
	using Vertex2 = GMap3::Vertex2;
	using HalfEdge = GMap3::HalfEdge;
	using Edge = GMap3::Edge;
	using Edge2 = GMap3::Edge2;
	using Face = GMap3::Face;
	using Face2 = GMap3::Face2;
	using Volume = GMap3::Volume;

	using Cells = std::tuple<Vertex, Vertex2, HalfEdge, Edge, Edge2, Face, Face2, Volume>;
	static constexpr const char* cell_names[] = {"Vertex", "Vertex2", "HalfEdge", "Edge",
												 "Edge2",  "Face",	  "Face2",	  "Volume"};

	template <typename T>
	using Attribute = MapBase::Attribute<T>;
	using AttributeGen = MapBase::AttributeGen;
	using MarkAttribute = MapBase::MarkAttribute;
};


GMap3::Vertex CGOGN_CORE_EXPORT cut_edge(GMap3& m, GMap3::Edge e, bool set_indices = true);

GMap3::Edge cut_face(GMap3& m, GMap3::Vertex v1, GMap3::Vertex v2, bool set_indices = true);

GMap3::Face cut_volume(GMap3& m, const std::vector<Dart>& path, bool set_indices = true);

GMap3::Volume close_hole(GMap3& m, Dart d, bool set_indices);

uint32 close(GMap3& m, bool set_indices);

inline Dart beta3(const GMap3& m, Dart d)
{
	return (*(m.beta3_))[d.index];
}


inline void beta3_sew(GMap3& m, Dart d, Dart e)
{
	cgogn_assert(beta3(m, d) == d);
	cgogn_assert(beta3(m, e) == e);
	(*(m.beta3_))[d.index] = e;
	(*(m.beta3_))[e.index] = d;
}

inline void beta3_unsew(GMap3& m, Dart d)
{
	Dart e = beta3(m, d);
	(*(m.beta3_))[d.index] = d;
	(*(m.beta3_))[e.index] = e;
}

inline Dart phi3(const GMap3& m, Dart d)
{
	Dart dd = beta3(m, d);
	if (dd == d)
		return d;
	return beta0(m, dd);
}

inline void phi3_sew(GMap3& m, Dart d, Dart e)
{
	beta3_sew(m, d, beta0(m, e));
	beta3_sew(m, beta0(m, d), e);
}

inline void phi3_unsew(GMap3& m, Dart d)
{
	Dart e = beta0(m, d);
	beta3_unsew(m, d);
	beta3_unsew(m, e);
}

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_GMAP_GMAP3_H_

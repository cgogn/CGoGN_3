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

#ifndef CGOGN_CORE_TYPES_GMAP_GMAP2_H_
#define CGOGN_CORE_TYPES_GMAP_GMAP2_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/cmap/gmap/gmap1.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT GMap2 : public GMap1
{
	static const uint8 dimension = 2;

	using Vertex = Cell<Orbit::BETA1_BETA2>;
	using HalfEdge = Cell<Orbit::BETA0>;
	using Edge = Cell<Orbit::BETA0_BETA2>;
	using Face = Cell<Orbit::BETA0_BETA1>;
	using Volume = Cell<Orbit::BETA0_BETA1_BETA2>;
	using CC = Volume;

	using Cells = std::tuple<Vertex, HalfEdge, Edge, Face, Volume>;

	std::shared_ptr<Attribute<Dart>> beta2_;

	GMap2() : GMap1()
	{
		beta2_ = add_relation("beta2");
	}

};

template <>
struct mesh_traits<GMap2>
{
	using MeshType = GMap2;
	using ParentDimType = GMap1;
	static constexpr const char* name = "GMap2";
	static constexpr const uint8 dimension = 2;

	using Vertex = GMap2::Vertex;
	using HalfEdge = GMap2::HalfEdge;
	using Edge = GMap2::Edge;
	using Face = GMap2::Face;
	using Volume = GMap2::Volume;

	using Cells = std::tuple<Vertex, HalfEdge, Edge, Face, Volume>;
	static constexpr const char* cell_names[] = {"Vertex", "HalfEdge", "Edge", "Face", "Volume"};

	template <typename T>
	using Attribute = CMapBase::Attribute<T>;
	using AttributeGen = CMapBase::AttributeGen;
	using MarkAttribute = CMapBase::MarkAttribute;
};

GMap2::Vertex CGOGN_CORE_EXPORT cut_edge(GMap2& m, GMap2::Edge e, bool set_indices = true);

GMap2::Edge CGOGN_CORE_EXPORT cut_face(GMap2& m, GMap2::Vertex v1, GMap2::Vertex v2, bool set_indices = true);

bool CGOGN_CORE_EXPORT flip_edge(GMap2& m, GMap2::Edge e, bool set_indices = true);

void CGOGN_CORE_EXPORT merge_incident_faces(GMap2& m, GMap2::Edge e, bool set_indices = true);

GMap2::Face CGOGN_CORE_EXPORT add_face(GMap2& m, uint32 size, bool set_indices = true);

GMap2::Volume CGOGN_CORE_EXPORT add_pyramid(GMap2& m, uint32 size, bool set_indices = true);

GMap2::Volume CGOGN_CORE_EXPORT add_prism(GMap2& m, uint32 size, bool set_indices = true);

void CGOGN_CORE_EXPORT remove_volume(GMap2& m, GMap2::Volume v);

//void CGOGN_CORE_EXPORT reverse_orientation(GMap2& m);

bool CGOGN_CORE_EXPORT edge_can_collapse(const GMap2& m, GMap2::Edge e);

bool CGOGN_CORE_EXPORT edge_can_flip(const GMap2& m, GMap2::Edge e);

bool CGOGN_CORE_EXPORT check_integrity(GMap2& m, bool verbose = true);

GMap2::Face close_hole(GMap2& m, Dart d, bool set_indices);

int32 close(GMap2& m, bool set_indices=false);


inline Dart beta2(const GMap2& m, Dart d)
{
	return (*(m.beta2_))[d.index];
}

inline Dart boundary_beta1(const GMap2& m, Dart d)
{
	return beta2(m, d);
}

inline bool on_boundary(const GMap2& m, Dart d)
{
	return beta2(m, d) == d;
}


inline void beta2_sew(GMap2& m, Dart d, Dart e)
{
	cgogn_assert(beta2(m, d) == d);
	cgogn_assert(beta2(m, e) == e);
	(*(m.beta2_))[d.index] = e;
	(*(m.beta2_))[e.index] = d;
}

inline void beta2_unsew(GMap2& m, Dart d)
{
	Dart e = beta2(m, d);
	(*(m.beta2_))[d.index] = d;
	(*(m.beta2_))[e.index] = e;
}

inline Dart phi2(const GMap2& m, Dart d)
{
	Dart dd = beta2(m, d);
	if (dd == d)
		return d;
	return beta0(m, dd);
}

inline void phi2_sew(GMap2& m, Dart d, Dart e)
{
	beta2_sew(m, d, beta0(m,e));
	beta2_sew(m, beta0(m, d), e);
}

inline void phi2_unsew(GMap2& m, Dart d)
{
	Dart e = beta0(m, d);
	beta2_unsew(m, d);
	beta2_unsew(m, e);
}


} // namespace cgogn

#endif // CGOGN_CORE_TYPES_GMAP_GMAP2_H_

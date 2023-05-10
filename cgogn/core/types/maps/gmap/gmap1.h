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

#ifndef CGOGN_CORE_TYPES_GMAP_GMAP1_H_
#define CGOGN_CORE_TYPES_GMAP_GMAP1_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/maps/gmap/gmap0.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT GMap1 : public GMap0
{
	static const uint8 dimension = 1;

	using Parent = GMap0;

	using Vertex = Cell<Orbit::BETA1>;
	using Edge = typename GMap0::Edge;
	using Face = Cell<Orbit::BETA0_BETA1>;

	using Cells = std::tuple<Vertex, Edge, Face>;

	std::shared_ptr<Attribute<Dart>> beta1_;

	GMap1() : GMap0()
	{
		beta1_ = add_relation("beta1");
	}

};

template <>
struct mesh_traits<GMap1>
{
	using MeshType = GMap1;
	using Parent = GMap1::Parent;
	static constexpr const uint8 dimension = 1;

	using Vertex = GMap1::Vertex;
	using Edge = GMap1::Edge;
	using Face = GMap1::Face;

	using Cells = std::tuple<Vertex, Edge, Face>;
	static constexpr const char* cell_names[] = {"Vertex", "Edge", "Face"};

	template <typename T>
	using Attribute = MapBase::Attribute<T>;
	using AttributeGen = MapBase::AttributeGen;
	static constexpr const char* name = "GMap1";
	using MarkAttribute = MapBase::MarkAttribute;
};


GMap1::Vertex CGOGN_CORE_EXPORT cut_edge(GMap1& m, GMap1::Edge e, bool set_indices = true);

GMap1::Vertex CGOGN_CORE_EXPORT collapse_edge(GMap1& m, GMap1::Edge e, bool set_indices);

GMap1::Face CGOGN_CORE_EXPORT add_face(GMap1& m, uint32 size, bool set_indices = true);

void CGOGN_CORE_EXPORT remove_face(GMap1& m, GMap1::Face f);

bool CGOGN_CORE_EXPORT check_integrity(GMap1& m, bool verbose = true);


inline Dart beta1(const GMap1& m, Dart d)
{
	return (*(m.beta1_))[d.index];
}

inline Dart boundary_beta1(const GMap1& m, Dart d)
{
	return beta1(m,d);
}


inline Dart phi1(const GMap1& m, Dart d)
{
	return beta1(m,beta0(m,d));
}


inline Dart phi_1(const GMap1& m, Dart d)
{
	return beta0(m,beta1(m,d));
}


inline void beta1_sew(GMap1& m, Dart d, Dart e)
{
	cgogn_assert(beta1(m, d) == d);
	cgogn_assert(beta1(m, e) == e);
	(*(m.beta1_))[d.index] = e;
	(*(m.beta1_))[e.index] = d;
}

inline void beta1_unsew(GMap1& m, Dart d)
{
	Dart e = beta1(m, d);
	(*(m.beta1_))[d.index] = d;
	(*(m.beta1_))[e.index] = e;
}




/*}****************************************************************************/

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_GMAP_GMAP1_H_


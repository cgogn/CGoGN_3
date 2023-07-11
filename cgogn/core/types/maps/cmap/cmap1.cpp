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

#include <cgogn/core/types/maps/cmap/cmap1.h>

namespace cgogn
{

CMap1::Vertex cut_edge(CMap1& m, CMap1::Edge e, bool set_indices)
{
	Dart d = add_dart(m);
	phi1_sew(m, e.dart, d);
	CMap1::Vertex v(d);

	if (set_indices)
	{
		if (is_indexed<CMap1::Vertex>(m))
			set_index(m, v, new_index<CMap1::Vertex>(m));
		// CMap1::Edge is the same orbit as CMap1::Vertex
		if (is_indexed<CMap1::Face>(m))
			copy_index<CMap1::Face>(m, d, e.dart);
	}

	return v;
}

CMap1::Vertex collapse_edge(CMap1& m, CMap1::Edge e, bool set_indices)
{
	Dart d = phi_1(m, e.dart);
	phi1_unsew(m, d);
	remove_dart(m, e.dart);
	CMap1::Vertex v(phi1(m, d));

	if (set_indices)
	{
	}

	return v;
}

CMap1::Face add_face(CMap1& m, uint32 size, bool set_indices)
{
	Dart d = add_dart(m);
	for (uint32 i = 1u; i < size; ++i)
	{
		Dart e = add_dart(m);
		phi1_sew(m, d, e);
	}
	CMap1::Face f(d);

	if (set_indices)
	{
		if (is_indexed<CMap1::Vertex>(m))
		{
			foreach_incident_vertex(
				m, f,
				[&](CMap1::Vertex v) -> bool {
					set_index(m, v, new_index<CMap1::Vertex>(m));
					return true;
				},
				MapBase::TraversalPolicy::DART_MARKING);
		}
		// CMap1::Edge is the same orbit as CMap1::Vertex
		if (is_indexed<CMap1::Face>(m))
			set_index(m, f, new_index<CMap1::Face>(m));
	}

	return f;
}

void remove_face(CMap1& m, CMap1::Face f)
{
	Dart it = phi1(m, f.dart);
	while (it != f.dart)
	{
		Dart next = phi1(m, it);
		remove_dart(m, it);
		it = next;
	}
	remove_dart(m, f.dart);
}

bool check_integrity(CMap1& m, bool verbose)
{
	bool result = true;
	for (Dart d = m.begin(), end = m.end(); d != end; d = m.next(d))
	{
		bool relations = phi<-1, 1>(m, d) == d && phi<1, -1>(m, d) == d;
		if (verbose && !relations)
			std::cerr << "Dart " << d << " has bad relations" << std::endl;

		result &= relations;
	}
	result &= check_indexing<CMap1::Vertex>(m);
	result &= check_indexing<CMap1::Edge>(m);
	result &= check_indexing<CMap1::Face>(m);
	return result;
}

void phi1_sew(CMap1& m, Dart d, Dart e)
{
	Dart f = phi1(m, d);
	Dart g = phi1(m, e);
	(*(m.phi1_))[d.index] = g;
	(*(m.phi1_))[e.index] = f;
	(*(m.phi_1_))[g.index] = d;
	(*(m.phi_1_))[f.index] = e;
}

void phi1_unsew(CMap1& m, Dart d)
{
	Dart e = phi1(m, d);
	Dart f = phi1(m, e);
	(*(m.phi1_))[d.index] = f;
	(*(m.phi1_))[e.index] = e;
	(*(m.phi_1_))[f.index] = d;
	(*(m.phi_1_))[e.index] = e;
}

} // namespace cgogn

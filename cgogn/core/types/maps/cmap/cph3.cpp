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

#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/types/maps/cmap/cph3.h>

#include <cgogn/core/types/cell_marker.h>

#include <unordered_set>

namespace cgogn
{

/***************************************************
 *              LEVELS MANAGEMENT                  *
 ***************************************************/

uint32 CPH3::dart_level(Dart d) const
{
	return (*dart_level_)[d.index];
}

void CPH3::set_dart_level(Dart d, uint32 l)
{
	if (uint32(nb_darts_per_level_.size()) > dart_level(d))
		nb_darts_per_level_[dart_level(d)]--;
	if (uint32(nb_darts_per_level_.size()) < l)
		nb_darts_per_level_.resize(l);
	nb_darts_per_level_[l]++;
	if (l > dart_level(d) && l > maximum_level_)
		maximum_level_ = l;
	if (l < dart_level(d))
	{
		while (nb_darts_per_level_[maximum_level_] == 0u)
		{
			--maximum_level_;
			nb_darts_per_level_.pop_back();
		}
	}
	(*dart_level_)[d.index] = l;
}

/***************************************************
 *             EDGE ID MANAGEMENT                  *
 ***************************************************/

uint32 CPH3::edge_id(Dart d) const
{
	return (*edge_id_)[d.index];
}

void CPH3::set_edge_id(Dart d, uint32 i)
{
	(*edge_id_)[d.index] = i;
}

uint32 CPH3::refinement_edge_id(Dart d, Dart e) const
{
	uint32 d_id = edge_id(d);
	uint32 e_id = edge_id(e);

	uint32 id = d_id + e_id;

	if (id == 0u)
		return 1u;
	else if (id == 1u)
		return 2u;
	else if (id == 2u)
	{
		if (d_id == e_id)
			return 0u;
		else
			return 1u;
	}
	// else if (id == 3)
	return 0u;
}

/***************************************************
 *             FACE ID MANAGEMENT                  *
 ***************************************************/

uint32 CPH3::face_id(Dart d) const
{
	return (*face_id_)[d.index];
}

void CPH3::set_face_id(Dart d, uint32 i)
{
	(*face_id_)[d.index] = i;
}

uint32 CPH3::refinement_face_id(const std::vector<Dart>& cut_path) const
{
	std::unordered_set<uint32> set_fid;
	for (Dart d : cut_path)
		set_fid.insert(face_id(d));
	uint32 result = 0;
	while (set_fid.find(result) != set_fid.end())
		++result;
	return result;
}

/***************************************************
 *                  EDGE INFO                      *
 ***************************************************/

uint32 CPH3::edge_level(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	return std::max(dart_level(d), dart_level(phi1(*this, d)));
}

Dart CPH3::edge_youngest_dart(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	Dart d2 = phi2(*this, d);
	if (dart_level(d) > dart_level(d2))
		return d;
	return d2;
}

bool CPH3::edge_is_subdivided(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	if (current_level_ == maximum_level_)
		return false;

	Dart d1 = phi1(*this, d);
	CPH3 m(*this);
	m.current_level_++;
	Dart d1_l = phi1(m, d);
	if (d1 != d1_l)
		return true;
	else
		return false;
}

/***************************************************
 *                  FACE INFO                      *
 ***************************************************/

uint32 CPH3::face_level(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	if (current_level_ == 0)
		return 0;

	Dart it = d;
	Dart old = it;
	uint32 l_old = dart_level(old);
	uint32 fLevel = edge_level(it);
	do
	{
		it = phi1(*this, it);
		uint32 dl = dart_level(it);

		// compute the oldest dart of the face in the same time
		if (dl < l_old)
		{
			old = it;
			l_old = dl;
		}
		uint32 l = edge_level(it);
		fLevel = l < fLevel ? l : fLevel;
	} while (it != d);

	CPH3 m(*this);
	m.current_level_ = fLevel;

	uint32 nbSubd = 0;
	it = old;
	uint32 eId = m.edge_id(old);
	uint32 init_dart_level = m.dart_level(it);
	do
	{
		++nbSubd;
		it = phi1(m, it);
	} while (m.edge_id(it) == eId && m.dart_level(it) != init_dart_level);

	while (nbSubd > 1)
	{
		nbSubd /= 2;
		--fLevel;
	}

	return fLevel;
}

Dart CPH3::face_origin(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	Dart p = d;
	uint32 pLevel = dart_level(p);
	CPH3 m(*this);
	do
	{
		m.current_level_ = pLevel;
		p = m.face_oldest_dart(p);
		pLevel = m.dart_level(p);
	} while (pLevel > 0);
	return p;
}

Dart CPH3::face_oldest_dart(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	Dart it = d;
	Dart oldest = it;
	uint32 l_old = dart_level(oldest);
	do
	{
		uint32 l = dart_level(it);
		if (l == 0)
			return it;

		if (l < l_old)
		//		if(l < l_old || (l == l_old && it < oldest))
		{
			oldest = it;
			l_old = l;
		}
		it = phi1(*this, it);
	} while (it != d);

	return oldest;
}

Dart CPH3::face_youngest_dart(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	Dart it = d;
	Dart youngest = it;
	uint32 l_young = dart_level(youngest);
	do
	{
		uint32 l = dart_level(it);
		if (l == current_level_)
			return it;

		if (l > l_young)
		//		if(l < l_young || (l == l_young && it < youngest))
		{
			youngest = it;
			l_young = l;
		}
		it = phi1(*this, it);
	} while (it != d);

	return youngest;
}

bool CPH3::face_is_subdivided(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");
	uint32 fLevel = face_level(d);
	if (fLevel < current_level_)
		return false;

	bool subd = false;
	CPH3 m(*this);
	m.current_level_++;
	if (m.dart_level(phi1(m, d)) == current_level_ && m.edge_id(phi1(m, d)) != m.edge_id(d))
		subd = true;
	return subd;
}

bool CPH3::face_is_subdivided_once(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");

	uint32 fLevel = face_level(d);
	if (fLevel < current_level_)
		return false;

	uint32 degree = 0;
	bool subd = false;
	bool subdOnce = true;
	Dart fit = d;
	CPH3 m(*this), m2(*this);
	m.current_level_ = current_level_ + 1;
	m2.current_level_ = current_level_ + 2;
	do
	{
		if (dart_level(phi1(m, fit)) == current_level_ && edge_id(phi1(m, fit)) != edge_id(fit))
		{
			subd = true;
			if (dart_level(phi1(m2, fit)) == current_level_ && edge_id(phi1(m2, fit)) != edge_id(fit))
				subdOnce = false;
		}
		++degree;
		fit = phi1(*this, fit);
	} while (subd && subdOnce && fit != d);

	if (degree == 3 && subd)
	{
		Dart cf = phi2(m, phi1(m, d));
		if (dart_level(phi1(m2, cf)) == current_level_ && edge_id(phi1(m2, cf)) != edge_id(cf))
			subdOnce = false;
	}

	return subd && subdOnce;
}

/***************************************************
 *                 VOLUME INFO                     *
 ***************************************************/

uint32 CPH3::volume_level(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");

	// The level of a volume is the
	// minimum of the levels of its faces

	Dart oldest = d;
	uint32 lold = dart_level(oldest);
	uint32 vLevel = std::numeric_limits<uint32>::max();

	foreach_incident_face(*this, CPH3::CMAP::Volume(d), [&](CPH3::CMAP::Face f) -> bool {
		uint32 fLevel = face_level(f.dart);
		vLevel = fLevel < vLevel ? fLevel : vLevel;
		Dart old = face_oldest_dart(f.dart);
		if (dart_level(old) < lold)
		{
			oldest = old;
			lold = dart_level(old);
		}
		return true;
	});

	CPH3 m(*this);
	m.current_level_ = vLevel;

	uint32 nbSubd = 0;
	Dart it = oldest;
	uint32 eId = edge_id(oldest);
	do
	{
		++nbSubd;
		it = phi<1, 2, 1>(m, it);
	} while (edge_id(it) == eId && lold != dart_level(it));

	while (nbSubd > 1)
	{
		nbSubd /= 2;
		--vLevel;
	}

	return vLevel;
}

Dart CPH3::volume_oldest_dart(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");

	Dart oldest = d;
	uint32 l_old = dart_level(oldest);
	foreach_incident_face(*this, CPH3::CMAP::Volume(oldest), [&](CPH3::CMAP::Face f) -> bool {
		Dart old = face_oldest_dart(f.dart);
		uint32 l = dart_level(old);
		if (l < l_old)
		{
			oldest = old;
			l_old = l;
		}
		return true;
	});

	return oldest;
}

Dart CPH3::volume_youngest_dart(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");

	Dart youngest = d;
	uint32 l_young = dart_level(youngest);
	foreach_incident_face(*this, CPH3::CMAP::Volume(youngest), [&](CPH3::CMAP::Face f) -> bool {
		Dart young = face_youngest_dart(f.dart);
		uint32 l = dart_level(young);
		if (l > l_young)
		{
			youngest = young;
			l_young = l;
		}
		return true;
	});

	return youngest;
}

bool CPH3::volume_is_subdivided(Dart d) const
{
	cgogn_message_assert(dart_level(d) <= current_level_, "Access to a dart introduced after current level");

	uint32 vLevel = volume_level(d);
	if (vLevel < current_level_)
		return false;

	bool faceAreSubdivided = face_is_subdivided(d);

	foreach_incident_face(*this, CPH3::CMAP::Volume(d), [&](CPH3::CMAP::Face f) -> bool {
		faceAreSubdivided &= face_is_subdivided(f.dart);
		return true;
	});

	bool subd = false;
	CPH3 m(*this);
	m.current_level_++;
	if (faceAreSubdivided && dart_level(phi<1, 1, 2>(m, d)) == current_level_ &&
		face_id(phi<1, 1, 2>(m, d)) != face_id(d))
		subd = true;

	return subd;
}

/*************************************************************************/
// Basic phi functions
/*************************************************************************/

Dart phi2bis(const CPH3& m, Dart d)
{
	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);

	uint32 face_id = m.face_id(d);
	Dart it = d;

	it = phi2(map, it);

	if (m.face_id(it) == face_id)
		return it;
	else
	{
		do
		{
			it = phi2(map, phi3(map, it));
		} while (m.face_id(it) != face_id);

		return it;
	}
}

Dart phi1(const CPH3& m, Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);

	if (m.current_level_ == m.maximum_level_)
		return phi1(map, d);

	bool finished = false;
	uint32 edge_id = m.edge_id(d);
	Dart it = d;
	do
	{
		it = phi1(map, it);
		if (m.dart_level(it) <= m.current_level_)
			finished = true;
		else
			while (m.edge_id(it) != edge_id)
				it = phi1(map, phi2bis(m, it));
	} while (!finished);
	return it;
}

Dart phi_1(const CPH3& m, Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);

	bool finished = false;
	Dart it = phi_1(map, d);
	uint32 edge_id = m.edge_id(it);

	do
	{
		if (m.dart_level(it) <= m.current_level_)
			finished = true;
		else
		{
			it = phi_1(map, it);
			while (m.edge_id(it) != edge_id)
				it = phi_1(map, phi2bis(m, it));
		}
	} while (!finished);
	return it;
}

Dart phi2(const CPH3& m, Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);
	return phi2(map, phi_1(map, phi1(m, d)));
}

Dart phi3(const CPH3& m, Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);
	if (phi3(map, d) == d)
		return d;
	return phi3(map, phi_1(map, phi1(m, d)));
}

/*************************************************************************/
// Specific overloads
/*************************************************************************/

Dart add_dart(CPH3& m)
{
	Dart d = add_dart(static_cast<CPH3::CMAP&>(m));
	if (uint32(m.nb_darts_per_level_.size()) < m.current_level_)
		m.nb_darts_per_level_.resize(m.current_level_);
	m.nb_darts_per_level_[m.current_level_]++;
	m.set_edge_id(d, 0u);
	m.set_face_id(d, 0u);
	m.set_dart_level(d, m.current_level_);

	// update max level if needed
	if (m.current_level_ > m.maximum_level_)
		m.maximum_level_ = m.current_level_;

	return d;
}

/*************************************************************************/
// Operators
/*************************************************************************/

CPH3::CMAP::Vertex cut_edge(CPH3& m, CPH3::CMAP::Edge e, bool set_indices)
{
	CPH3::CMAP& map = static_cast<CPH3::CMAP&>(m);

	CPH3::CMAP::Vertex v = cut_edge(map, e, false);

	Dart d = e.dart;
	do
	{
		m.set_edge_id(phi1(map, d), m.edge_id(d));
		m.set_edge_id(phi3(map, d), m.edge_id(d));
		m.set_edge_id(phi2(map, d), m.edge_id(phi<1, 2>(map, d)));
		m.set_face_id(phi1(map, d), m.face_id(d));
		m.set_face_id(phi3(map, d), m.face_id(d));
		m.set_face_id(phi2(map, d), m.face_id(phi<1, 2>(map, d)));
		m.set_dart_level(phi1(map, d), m.current_level_);
		m.set_dart_level(phi2(map, d), m.current_level_);
		d = phi<2, 3>(map, d);
	} while (d != e.dart);
	if (set_indices)
	{
		if (is_indexed<CPH3::CMAP::Vertex>(m))
			set_index(map, v, new_index<CPH3::CMAP::Vertex>(m));
		if (is_indexed<CPH3::CMAP::Edge>(m))
		{
			uint32 ne = new_index<CPH3::CMAP::Edge>(m);
			foreach_dart_of_orbit(map, e, [&](Dart d) -> bool {
				if (m.dart_level(d) == m.current_level_)
					set_index<CPH3::CMAP::Edge>(m, d, ne);
				return true;
			});
			ne = new_index<CPH3::CMAP::Edge>(m);
			foreach_dart_of_orbit(map, CPH3::CMAP::Edge(phi1(map, e.dart)), [&](Dart d) -> bool {
				if (m.dart_level(d) == m.current_level_)
					set_index<CPH3::CMAP::Edge>(m, d, ne);
				return true;
			});
		}
		if (is_indexed<CPH3::CMAP::Face>(m))
		{
			d = e.dart;
			do
			{
				Dart it = phi1(map, d);
				do
				{
					it = phi1(map, it);
				} while (m.dart_level(it) < m.current_level_ - 1 && it != d);

				copy_index<CPH3::CMAP::Face>(map, phi1(map, d), it);
				it = phi2(map, d);
				do
				{
					it = phi1(map, it);
				} while (m.dart_level(it) < m.current_level_ - 1 && it != phi2(map, phi1(map, d)));
				copy_index<CPH3::CMAP::Face>(map, phi2(map, d), it);
				d = phi<2, 3>(map, d);
			} while (d != e.dart);
		}
		if (is_indexed<CPH3::CMAP::Volume>(m))
		{
			d = e.dart;
			do
			{
				if (!is_boundary(m, d))
				{
					Dart it = phi1(map, d);
					do
					{
						it = phi1(map, it);
					} while (m.dart_level(it) < m.current_level_ - 1 && it != d);
					copy_index<CPH3::CMAP::Volume>(map, phi1(map, d), it);
					it = phi2(map, d);
					do
					{
						it = phi1(map, it);
					} while (m.dart_level(it) < m.current_level_ - 1 && it != phi2(map, phi1(map, d)));
					copy_index<CPH3::CMAP::Volume>(map, phi2(map, d), it);
				}
				d = phi<2, 3>(map, d);
			} while (d != e.dart);
		}
	}

	return v;
}

CPH3::CMAP::Edge cut_face(CPH3& m, CPH3::CMAP::Vertex v1, CPH3::CMAP::Vertex v2, bool set_indices)
{
	CPH3::CMAP& map = static_cast<CPH3::CMAP&>(m);

	Dart d = v1.dart;
	Dart e = v2.dart;

	Dart dd = phi<3, 1>(m, v1.dart);
	Dart ee = phi<3, 1>(m, e);

	CPH3::CMAP::Edge result = cut_face(map, v1, v2, false);

	uint32 eid = m.refinement_edge_id(v1.dart, v2.dart);

	foreach_dart_of_orbit(m, result, [&](Dart d) -> bool {
		m.set_edge_id(d, eid);
		m.set_face_id(d, m.face_id(v1.dart));
		m.set_dart_level(d, m.current_level_);
		return true;
	});

	if (set_indices)
	{
		if (is_indexed<CPH3::CMAP::Vertex>(m))
		{
			copy_index<CPH3::CMAP::Vertex>(m, phi_1(m, e), v1.dart);
			copy_index<CPH3::CMAP::Vertex>(m, phi_1(m, ee), v1.dart);
			copy_index<CPH3::CMAP::Vertex>(m, phi_1(m, d), e);
			copy_index<CPH3::CMAP::Vertex>(m, phi_1(m, dd), e);
		}
		if (is_indexed<CPH3::CMAP::Edge>(m))
			set_index(m, CPH3::CMAP::Edge(phi_1(m, v1.dart)), new_index<CPH3::CMAP::Edge>(m));
		if (is_indexed<CPH3::CMAP::Face>(m))
		{
			uint32 nf1 = new_index<CPH3::CMAP::Face>(m);
			uint32 nf2 = new_index<CPH3::CMAP::Face>(m);
			foreach_dart_of_orbit(m, CPH3::CMAP::Face(d), [&](Dart df) -> bool {
				if (m.current_level_ == m.dart_level(df))
					set_index<CPH3::CMAP::Face>(m, df, nf1);
				return true;
			});
			foreach_dart_of_orbit(m, CPH3::CMAP::Face(e), [&](Dart df) -> bool {
				if (m.current_level_ == m.dart_level(df))
					set_index<CPH3::CMAP::Face>(m, df, nf2);
				return true;
			});
		}
		if (is_indexed<CPH3::CMAP::Volume>(m))
		{
			if (!is_boundary(m, d))
			{
				copy_index<CPH3::CMAP::Volume>(map, phi_1(m, d), d);
				copy_index<CPH3::CMAP::Volume>(map, phi_1(m, e), d);
			}
			if (!is_boundary(m, dd))
			{
				copy_index<CPH3::CMAP::Volume>(map, phi_1(m, dd), dd);
				copy_index<CPH3::CMAP::Volume>(map, phi_1(m, ee), dd);
			}
		}
	}

	return result;
}

CPH3::CMAP::Face cut_volume(CPH3& m, const std::vector<Dart>& path, bool set_indices)
{
	CPH3::CMAP& map = static_cast<CPH3::CMAP&>(m);

	uint32 vid = m.refinement_face_id(path);
	uint32 vlevel = m.volume_level(path[0]);

	CPH3::CMAP::Face result = cut_volume(map, path, false);

	Dart f0 = result.dart;
	Dart f1 = phi3(m, f0);

	foreach_dart_of_orbit(m, result, [&](Dart d) -> bool {
		m.set_edge_id(d, m.edge_id(phi2(m, d)));
		m.set_face_id(d, vid);
		m.set_dart_level(d, m.current_level_);
		return true;
	});

	if (set_indices)
	{
		if (is_indexed<CPH3::CMAP::Vertex>(m))
		{
			foreach_dart_of_orbit(m, CPH3::CMAP::Face(f0), [&](Dart d) -> bool {
				copy_index<CPH3::CMAP::Vertex>(m, d, phi<2, 1>(m, d));
				return true;
			});
		}
		if (is_indexed<CPH3::CMAP::Edge>(m))
		{
			foreach_dart_of_orbit(m, CPH3::CMAP::Face2(f0), [&](Dart d) -> bool {
				copy_index<CPH3::CMAP::Edge>(map, d, phi2(m, d));
				copy_index<CPH3::CMAP::Edge>(map, phi3(m, d), d);
				return true;
			});
		}
		if (is_indexed<CPH3::CMAP::Face>(m))
			set_index(m, CPH3::CMAP::Face(f0), new_index<CPH3::CMAP::Face>(m));
		if (is_indexed<CPH3::CMAP::Volume>(m))
		{
			if (vlevel == m.current_level_)
			{
				foreach_dart_of_orbit(m, CPH3::CMAP::Face2(f0), [&](Dart d) -> bool {
					copy_index<CPH3::CMAP::Volume>(map, d, phi2(m, d));
					return true;
				});
				set_index(m, CPH3::CMAP::Volume(f1), new_index<CPH3::CMAP::Volume>(m));
			}
			else
			{
				uint32 ved1 = new_index<CPH3::CMAP::Volume>(m);
				uint32 ved2 = new_index<CPH3::CMAP::Volume>(m);
				foreach_dart_of_orbit(m, CPH3::CMAP::Volume(f0), [&](Dart d) -> bool {
					if (m.dart_level(d) == m.current_level_)
						set_index<CPH3::CMAP::Volume>(m, d, ved1);
					return true;
				});
				foreach_dart_of_orbit(m, CPH3::CMAP::Volume(f1), [&](Dart d) -> bool {
					if (m.dart_level(d) == m.current_level_)
						set_index<CPH3::CMAP::Volume>(m, d, ved2);
					return true;
				});
			}
		}
	}

	return result;
}

} // namespace cgogn

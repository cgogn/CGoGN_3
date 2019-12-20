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

#ifndef CGOGN_CORE_TYPES_CMAP_CPH3_H_
#define CGOGN_CORE_TYPES_CMAP_CPH3_H_

#include <cgogn/core/types/cmap/cmap3.h>

#include <unordered_set>

namespace cgogn
{

struct CPH3
{
	using CMAP = CMap3;

	template <typename T>
	using Attribute = CMAP::Attribute<T>;
	using AttributeGen = CMAP::AttributeGen;
	using MarkAttribute = CMAP::MarkAttribute;

	CMAP& m_;

	std::shared_ptr<Attribute<uint32>> dart_level_;
	std::shared_ptr<Attribute<uint32>> edge_id_;
	std::shared_ptr<Attribute<uint32>> face_id_;

	std::vector<uint32>& nb_darts_per_level_;
	uint32& maximum_level_;

	uint32 current_level_;

	CPH3(CMAP& m)
		: m_(m), current_level_(0),
		  nb_darts_per_level_(m.get_attribute<std::vector<uint32>>("cph3_nb_darts_per_level")),
		  maximum_level_(m.get_attribute<uint32>("cph3_maximum_level"))
	{
		dart_level_ = m_.darts_.get_attribute<uint32>("dart_level");
		if (!dart_level_)
			dart_level_ = m_.darts_.add_attribute<uint32>("dart_level");
		else
		{
			for (uint32 l : *dart_level_)
				if (l > maximum_level_)
					maximum_level_ = l;
		}

		edge_id_ = m_.darts_.get_attribute<uint32>("edge_id");
		if (!edge_id_)
			edge_id_ = m_.darts_.add_attribute<uint32>("edge_id");

		face_id_ = m_.darts_.get_attribute<uint32>("face_id");
		if (!face_id_)
			face_id_ = m_.darts_.add_attribute<uint32>("face_id");
	}

	CPH3(const CPH3& cph3)
		: m_(cph3.m_), dart_level_(cph3.dart_level_), edge_id_(cph3.edge_id_), face_id_(cph3.face_id_),
		  nb_darts_per_level_(cph3.nb_darts_per_level_), maximum_level_(cph3.maximum_level_),
		  current_level_(cph3.current_level_)
	{
	}

	operator CMAP&()
	{
		return m_;
	}
	operator const CMAP&() const
	{
		return m_;
	}

	inline Dart begin() const
	{
		Dart d(m_.darts_.first_index());
		uint32 lastidx = m_.darts_.last_index();
		while (dart_level(d) > current_level_ && d.index < lastidx)
			d = Dart(m_.darts_.next_index(d.index));
		return d;
	}

	inline Dart end() const
	{
		return Dart(m_.darts_.last_index());
	}

	inline Dart next(Dart d) const
	{
		uint32 lastidx = m_.darts_.last_index();
		do
		{
			d = Dart(m_.darts_.next_index(d.index));
		} while (dart_level(d) > current_level_ && d.index < lastidx);
		return d;
	}

	/***************************************************
	 *              LEVELS MANAGEMENT                  *
	 ***************************************************/

	uint32 dart_level(Dart d) const;
	void set_dart_level(Dart d, uint32 l);
	void change_dart_level(Dart d, uint32 l);

	/***************************************************
	 *             EDGE ID MANAGEMENT                  *
	 ***************************************************/

	uint32 edge_id(Dart d) const;
	void set_edge_id(Dart d, uint32 i);
	uint32 refinement_edge_id(Dart d, Dart e) const;

	/***************************************************
	 *             FACE ID MANAGEMENT                  *
	 ***************************************************/

	uint32 face_id(Dart d) const;
	void set_face_id(Dart d, uint32 i);
	uint32 refinement_face_id(const std::vector<Dart>& cut_path) const;

	/***************************************************
	 *                  EDGE INFO                      *
	 ***************************************************/

	uint32 edge_level(Dart d) const;
	Dart edge_youngest_dart(Dart d) const;

	/***************************************************
	 *                  FACE INFO                      *
	 ***************************************************/

	uint32 face_level(Dart d);
	Dart face_origin(Dart d);
	Dart face_oldest_dart(Dart d);
	Dart face_youngest_dart(Dart d);

	/***************************************************
	 *                 VOLUME INFO                     *
	 ***************************************************/

	uint32 volume_level(Dart d);
	Dart volume_oldest_dart(Dart d);
	Dart volume_youngest_dart(Dart d);
};

} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CMAP_CPH3_H_

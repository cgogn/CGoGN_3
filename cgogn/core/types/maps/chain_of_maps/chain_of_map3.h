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

#ifndef CGOGN_CORE_TYPES_CHAINOFMAP_CHMAP3_H_
#define CGOGN_CORE_TYPES_CHAINOFMAP_CHMAP3_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/types/maps/gmap/gmap2.h>
#include <cgogn/core/types/maps/chain_of_maps/chain_of_map2.h>

namespace cgogn
{
	
struct ChainOfMap3 : public ChainOfMap2
{
	static const uint8 dimension = 3;
	GMap2 gm_volumes_;
	using Volume = GMap2::Volume;

	std::shared_ptr<Attribute<Dart>> sigma3_0_;
	std::shared_ptr<Attribute<Dart>> sigma3_1_;
	std::shared_ptr<Attribute<Dart>> sigma3_2_;
	std::shared_ptr<Attribute<std::vector<Dart>>> sigma0_3_;
	std::shared_ptr<Attribute<std::vector<Dart>>> sigma1_3_;
	std::shared_ptr<Attribute<std::vector<Dart>>> sigma2_3_;

	ChainOfMap3() : ChainOfMap2()
	{
		sigma3_0_ = gm_volumes_.add_relation("sigma3_0");
		sigma3_1_ = gm_volumes_.add_relation("sigma3_1");
		sigma3_2_ = gm_volumes_.add_relation("sigma3_2");
		sigma0_3_ = relations_vect_.emplace_back(darts_.add_attribute<std::vector<Dart>>("sigma0_3_"));
		sigma1_3_ = relations_vect_.emplace_back(gm_edges_.darts_.add_attribute<std::vector<Dart>>("sigma1_2_"));
		sigma2_3_ = relations_vect_.emplace_back(gm_edges_.darts_.add_attribute<std::vector<Dart>>("sigma2_3_"));
	}
};

template <>
struct mesh_traits<ChainOfMap3>
{
	static constexpr const char* name = "ChMap3";
	static constexpr const uint8 dimension = 3;

	using Vertex = ChainOfMap0::Vertex;
	using Edge = ChainOfMap1::Edge;
	using Face = ChainOfMap2::Face;
	using Volume = ChainOfMap3::Volume;

	using Cells = std::tuple<Vertex, Edge, Face, Volume>;
	static constexpr const char* cell_names[] = {"Vertex", "Edge", "Face", "Volume"};

	template <typename T>
	using Attribute = ChainOfMapBase::Attribute<T>;
	using AttributeGen = ChainOfMapBase::AttributeGen;
	using MarkAttribute = ChainOfMapBase::MarkAttribute;
};

inline Dart sigma3_2(const ChainOfMap3& cm, Dart d)
{
	return (*(cm.sigma3_2_))[d.index];
}

const std::vector<Dart>& sigma2_3(const ChainOfMap3& cm, Dart d)
{
	return (*(cm.sigma2_3_))[d.index];
}



} // namespace cgogn

#endif // CGOGN_CORE_TYPES_CHAINOFMAP_CHMAP3_H_

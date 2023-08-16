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

#ifndef CGOGN_MODELING_DECIMATION_EDGE_QUEUE_SQEM_H_
#define CGOGN_MODELING_DECIMATION_EDGE_QUEUE_SQEM_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/types/slab_quadric.h>
#include <cgogn/geometry/types/vector_traits.h>
namespace cgogn
{

namespace modeling
{

template <typename NONMANIFOLD>
struct DecimationSQEM_Helper
{
	template <typename T>
	using Attribute = typename mesh_traits<NONMANIFOLD>::template Attribute<T>;

	using Vertex = typename NONMANIFOLD::Vertex;
	using Edge = typename NONMANIFOLD::Edge;
	using Face = typename NONMANIFOLD::Face;

	using Vec3 = typename geometry::Vec3;
	using Vec4 = typename geometry::Vec4;
	using Scalar = typename geometry::Scalar;
	using Slab_Quadric = typename geometry::Slab_Quadric;

	DecimationSQEM_Helper(NONMANIFOLD& m, const Attribute<Vec4>* sphere_info) : m_(m), sphere_info_(sphere_info)
	{
		sphere_quadric_ = add_attribute<Slab_Quadric, Vertex>(m, "__sphere_quadric");
		foreach_cell(m_, [&](Face f) -> bool {
			std::vector<Vertex> iv = incident_vertices(m_, f);
			Slab_Quadric q(value<Vec4>(m_, sphere_info_, iv[0]), value<Vec4>(m_, sphere_info_, iv[1]),
						   value<Vec4>(m_, sphere_info_, iv[2]));
			Slab_Quadric& q0 = value<Slab_Quadric>(m_, sphere_quadric_, iv[0]);
			Slab_Quadric& q1 = value<Slab_Quadric>(m_, sphere_quadric_, iv[1]);
			Slab_Quadric& q2 = value<Slab_Quadric>(m_, sphere_quadric_, iv[2]);
			value<Slab_Quadric>(m_, sphere_quadric_, iv[0]) += q;
			value<Slab_Quadric>(m_, sphere_quadric_, iv[1]) += q;
			value<Slab_Quadric>(m_, sphere_quadric_, iv[2]) += q;
			return true;
		});
	}
	~DecimationSQEM_Helper()
	{
		remove_attribute<Vertex>(m_, sphere_quadric_);
	}

	Scalar edge_cost(Edge e, const Vec4& p)
	{
		std::vector<Vertex> iv = incident_vertices(m_, e);
		Slab_Quadric q;
		Slab_Quadric& q0 = value<Slab_Quadric>(m_, sphere_quadric_, iv[0]);
		Slab_Quadric& q1 = value<Slab_Quadric>(m_, sphere_quadric_, iv[1]);
		q+= value<Slab_Quadric>(m_, sphere_quadric_, iv[0]);
		q+= value<Slab_Quadric>(m_, sphere_quadric_, iv[1]);
		return q.eval(p);
	}

	Vec4 edge_optimal(Edge e)
	{
		std::vector<Vertex> iv = incident_vertices(m_, e);
		Slab_Quadric q;
		q += value<Slab_Quadric>(m_, sphere_quadric_, iv[0]);
		q += value<Slab_Quadric>(m_, sphere_quadric_, iv[1]);
		Vec4 p;
		if (q.optimized(p))
			if (p[3] < 0)
				return Scalar(0.5) * (value<Vec4>(m_, sphere_info_, iv[0]) + value<Vec4>(m_, sphere_info_, iv[1]));
			else
				return p;
		else
			return Scalar(0.5) * (value<Vec4>(m_, sphere_info_, iv[0]) + value<Vec4>(m_, sphere_info_, iv[1]));
	}


	NONMANIFOLD& m_;
	const Attribute<Vec4>* sphere_info_;
	std::shared_ptr<Attribute<Slab_Quadric>> sphere_quadric_;
	Slab_Quadric q_;
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_DECIMATION_EDGE_QUEUE_SQEM_H_

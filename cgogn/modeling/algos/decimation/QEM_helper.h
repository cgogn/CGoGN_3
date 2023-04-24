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

#ifndef CGOGN_MODELING_DECIMATION_EDGE_QUEUE_QEM_H_
#define CGOGN_MODELING_DECIMATION_EDGE_QUEUE_QEM_H_


#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/utils/type_traits.h>
#include <cgogn/geometry/types/quadric.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{


namespace modeling
{

using geometry::Vec3;
using geometry::Scalar;
using geometry::Quadric;

template <typename MESH>
struct DecimationQEM_Helper
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	DecimationQEM_Helper(MESH& m, const Attribute<Vec3>* vertex_position) : m_(m), vertex_position_(vertex_position)
	{
		vertex_quadric_ = add_attribute<Quadric, Vertex>(m, "__vertex_quadric");
		parallel_foreach_cell(m_, [&](Vertex v) -> bool {
			value<Quadric>(m_, vertex_quadric_, v).zero();
			return true;
		});
		parallel_foreach_cell(m_, [&](Face f) -> bool {
			std::vector<Vertex> iv = incident_vertices(m_, f);
			Quadric q(value<Vec3>(m_, vertex_position_, iv[0]), value<Vec3>(m_, vertex_position_, iv[1]),
					  value<Vec3>(m_, vertex_position_, iv[2]));
			value<Quadric>(m_, vertex_quadric_, iv[0]) += q;
			value<Quadric>(m_, vertex_quadric_, iv[1]) += q;
			value<Quadric>(m_, vertex_quadric_, iv[2]) += q;
			return true;
		});
	}
	~DecimationQEM_Helper()
	{
		remove_attribute<Vertex>(m_, vertex_quadric_);
	}

	Scalar edge_cost(Edge e, const Vec3& p)
	{
		std::vector<Vertex> iv = incident_vertices(m_, e);
		Quadric q;
		q += value<Quadric>(m_, vertex_quadric_, iv[0]);
		q += value<Quadric>(m_, vertex_quadric_, iv[1]);
		return q.eval(p);
	}

	Vec3 edge_optimal(Edge e)
	{
		std::vector<Vertex> iv = incident_vertices(m_, e);
		Quadric q;
		q += value<Quadric>(m_, vertex_quadric_, iv[0]);
		q += value<Quadric>(m_, vertex_quadric_, iv[1]);
		Vec3 p;
		if (q.optimized(p))
			return p;
		else
			return Scalar(0.5) * (value<Vec3>(m_, vertex_position_, iv[0]) + value<Vec3>(m_, vertex_position_, iv[1]));
	}

	void before_collapse(Edge e)
	{
		std::vector<Vertex> iv = incident_vertices(m_, e);
		q_.zero();
		q_ += value<Quadric>(m_, vertex_quadric_, iv[0]);
		q_ += value<Quadric>(m_, vertex_quadric_, iv[1]);
	}

	void after_collapse(Vertex v)
	{
		value<Quadric>(m_, vertex_quadric_, v) = q_;
	}

	MESH& m_;
	const Attribute<Vec3>* vertex_position_;
	std::shared_ptr<Attribute<Quadric>> vertex_quadric_;
	Quadric q_;
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_DECIMATION_EDGE_QUEUE_QEM_H_

/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_MODELING_ALGOS_DECIMATION_EDGE_APPROXIMATOR_MID_EDGE_H_
#define CGOGN_MODELING_ALGOS_DECIMATION_EDGE_APPROXIMATOR_MID_EDGE_H_

#include <cgogn/modeling/algos/decimation/edge_approximator.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/edge.h>

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH, typename VEC>
class EdgeApproximator_MidEdge : public EdgeApproximator<MESH, VEC>
{
public:

	using Scalar = typename geometry::vector_traits<VEC>::Scalar;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	inline EdgeApproximator_MidEdge(
		const MESH& m,
		const typename mesh_traits<MESH>::template AttributePtr<VEC> position
	) : EdgeApproximator<MESH, VEC>(m),
		position_(position)
	{}
	virtual ~EdgeApproximator_MidEdge()
	{}

	void init()
	{}

	VEC operator()(Edge e) const
	{
		auto vertices = incident_vertices(this->m_, e);
		return Scalar(0.5) * (value<VEC>(this->m_, position_, vertices[0]) + value<VEC>(this->m_, position_, vertices[1]));
	}

private:

	const typename mesh_traits<MESH>::template AttributePtr<VEC> position_;
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_DECIMATION_EDGE_APPROXIMATOR_MID_EDGE_H_

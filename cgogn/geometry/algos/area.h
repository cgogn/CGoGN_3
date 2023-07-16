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

#ifndef CGOGN_GEOMETRY_ALGOS_AREA_H_
#define CGOGN_GEOMETRY_ALGOS_AREA_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/functions/area.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

// how to compute the local area of a vertex
enum VertexAreaPolicy
{
	THIRD,		// take the third of the area of neighboring faces
	BARYCENTER, // take the barycenter of each triangle
	VORONOI,	// take voronoi cell on each triangle
	MIXED		// take the voronoi cell clamped in the triangle
};

template <typename MESH>
Scalar area(const MESH& m, typename mesh_traits<MESH>::Face f,
			const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	std::vector<Vertex> vertices = incident_vertices(m, f);
	Scalar face_area{0};
	for (uint32 i = 1, size = uint32(vertices.size()); i < size - 1; ++i)
	{
		face_area += area(value<Vec3>(m, vertex_position, vertices[0]), value<Vec3>(m, vertex_position, vertices[i]),
						  value<Vec3>(m, vertex_position, vertices[(i + 1) % size]));
	}
	return face_area;
}

template <typename MESH>
Scalar area(const MESH& m, typename mesh_traits<MESH>::Vertex v,
			const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	return area(m, v, vertex_position, VertexAreaPolicy::BARYCENTER);
}

/**
 * @brief compute the local vertex area
 * @param[in] m mesh
 * @param[in] v the vertex to compute
 * @param[in] vertex_position position attribute of mesh
 * @param[in] area_policy the method to determine the local area of a vertex
 * @returns computed local area
 * @todo generalization for any polygon and neighboroud
 */
template <typename MESH>
Scalar area(const MESH& m, typename mesh_traits<MESH>::Vertex v,
			const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
			const VertexAreaPolicy area_policy)
{
	static_assert(mesh_traits<MESH>::dimension == 2, "MESH dimension should be 2");

	using Face = typename mesh_traits<MESH>::Face;
	using Vertex = typename mesh_traits<MESH>::Vertex;

	Scalar vertex_area{0};

	switch (area_policy)
	{
	case VertexAreaPolicy::BARYCENTER: {
		const Vec3& vertex = value<Vec3>(m, vertex_position, v);
		std::vector<Vertex> vertices = adjacent_vertices_through_edge(m, v);
		uint32 size = uint32(vertices.size());
		for (uint32 i = 0; i < size; ++i)
		{
			const Vec3& current_vertex = value<Vec3>(m, vertex_position, vertices[i]);
			const Vec3& next_vertex = value<Vec3>(m, vertex_position, vertices[(i + 1) % size]);

			Vec3 median = (current_vertex + vertex) * 0.5;
			Vec3 next_median = (next_vertex + vertex) * 0.5;
			Vec3 centroid = (median * 2. + next_vertex) / 3.; // median ratio

			vertex_area += area(median, centroid, vertex);
			vertex_area += area(next_median, centroid, vertex);
		}
	}
	case VertexAreaPolicy::VORONOI: {
		const Vec3& vertex = value<Vec3>(m, vertex_position, v);
		std::vector<Vertex> vertices = adjacent_vertices_through_edge(m, v);
		uint32 size = uint32(vertices.size());
		for (uint32 i = 0; i < size; ++i)
		{
			const Vec3& current_vertex = value<Vec3>(m, vertex_position, vertices[i]);
			const Vec3& next_vertex = value<Vec3>(m, vertex_position, vertices[(i + 1) % size]);

			Vec3 median = (current_vertex + vertex) * 0.5;
			Vec3 next_median = (next_vertex + vertex) * 0.5;
			Vec3 circumcenter = (current_vertex + next_vertex + vertex) / 3.;

			vertex_area += area(median, circumcenter, vertex);
			vertex_area += area(next_median, circumcenter, vertex);
		}
	}
	case VertexAreaPolicy::MIXED: {
		const Vec3& vertex = value<Vec3>(m, vertex_position, v);
		std::vector<Vertex> vertices = adjacent_vertices_through_edge(m, v);
		uint32 size = uint32(vertices.size());
		for (uint32 i = 0; i < size; ++i)
		{
			const Vec3& current_vertex = value<Vec3>(m, vertex_position, vertices[i]);
			const Vec3& next_vertex = value<Vec3>(m, vertex_position, vertices[(i + 1) % size]);

			Vec3 median = (current_vertex + vertex) * 0.5;
			Vec3 next_median = (next_vertex + vertex) * 0.5;
			Vec3 median_centroid = (current_vertex + next_vertex) * 0.5;		  // median between the 2 vertices
			Vec3 voronoi_centroid = (current_vertex + next_vertex + vertex) / 3.; // voronoi center
			Vec3 centroid = ((median_centroid - vertex).norm() < (voronoi_centroid - vertex).norm())
								? median_centroid
								: voronoi_centroid; // "clamp" the centroid in the triangle
			vertex_area += area(median, centroid, vertex);
			vertex_area += area(next_median, centroid, vertex);
		}
	}
	case VertexAreaPolicy::THIRD: {
		foreach_incident_face(m, v, [&](Face iface) -> bool {
			vertex_area += area(m, iface, vertex_position) / 3.0;
			return true;
		});
	}
	}

	return vertex_area;
}

/**
 * compute every vertex area
 */
template <typename CELL, typename MESH>
void compute_area(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
				  typename mesh_traits<MESH>::template Attribute<Scalar>* cell_area, VertexAreaPolicy area_policy)
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "MESH dimension should be >= 2");
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");

	using Vertex = typename mesh_traits<MESH>::Vertex;
	static_assert(std::is_same_v<CELL, Vertex>, "VertexAreaPolicy only supported with Vertex");

	parallel_foreach_cell(m, [&](Vertex c) -> bool {
		value<Scalar>(m, cell_area, c) = area(m, c, vertex_position, area_policy);
		return true;
	});
}

/**
 * compute every cell area
 */
template <typename CELL, typename MESH>
void compute_area(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
				  typename mesh_traits<MESH>::template Attribute<Scalar>* cell_area)
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "MESH dimension should be >= 2");
	static_assert(is_in_tuple_v<CELL, typename mesh_traits<MESH>::Cells>, "CELL not supported in this MESH");
	parallel_foreach_cell(m, [&](CELL c) -> bool {
		value<Scalar>(m, cell_area, c) = area(m, c, vertex_position);
		return true;
	});
}

/**
 * compute the mesh area
 */
template <typename MESH>
Scalar area(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Face = typename mesh_traits<MESH>::Face;

	Scalar result = 0;
	foreach_cell(m, [&](Face f) -> bool {
		result += area(m, f, vertex_position);
		return true;
	});
	return result;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_AREA_H_

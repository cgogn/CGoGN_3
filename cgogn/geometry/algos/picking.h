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

#ifndef CGOGN_GEOMETRY_ALGOS_PICKING_H_
#define CGOGN_GEOMETRY_ALGOS_PICKING_H_

#include <cgogn/core/functions/traversals/face.h>

#include <cgogn/geometry/functions/distance.h>
#include <cgogn/geometry/functions/intersection.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

namespace internal
{

template <typename MESH>
std::vector<std::tuple<typename mesh_traits<MESH>::Face, Vec3, Scalar>> picking(
	const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position, const Vec3& A,
	const Vec3& B)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;
	using SelectedFace = std::tuple<Face, Vec3, Scalar>;

	Vec3 AB = B - A;
	cgogn_message_assert(AB.squaredNorm() > 0.0, "line must be defined by 2 different points");
	AB.normalize();

	std::vector<std::vector<SelectedFace>> selected_per_thread(thread_pool()->nb_workers());
	// std::vector<std::vector<uint32>> ear_indices_per_thread(thread_pool()->nb_workers());

	parallel_foreach_cell(m, [&](Face f) -> bool {
		uint32 worker_index = current_worker_index();
		Vec3 intersection_point;
		std::vector<Vertex> vertices = incident_vertices(m, f);
		if (vertices.size() == 3)
		{
			if (intersection_ray_triangle(A, AB, value<Vec3>(m, vertex_position, vertices[0]),
										  value<Vec3>(m, vertex_position, vertices[1]),
										  value<Vec3>(m, vertex_position, vertices[2]), &intersection_point))
				selected_per_thread[worker_index].emplace_back(f, intersection_point,
															   (intersection_point - A).squaredNorm());
		}
		else
		{
			for (uint32 i = 0, size = uint32(vertices.size()); i + 2 < size; i++)
			{
				if (intersection_ray_triangle(A, AB, value<Vec3>(m, vertex_position, vertices[0]),
											  value<Vec3>(m, vertex_position, vertices[i + 1]),
											  value<Vec3>(m, vertex_position, vertices[i + 2]), &intersection_point))
				{
					selected_per_thread[worker_index].emplace_back(f, intersection_point,
																   (intersection_point - A).squaredNorm());
					break;
				}
			}
		}
		// else
		// {
		// 	std::vector<uint32>& ear_indices = ear_indices_per_thread[worker_index];
		// 	ear_indices.clear();
		// 	append_ear_triangulation(m, f, position, ear_indices);
		// 	for (std::size_t i = 0; i < uint32(ear_indices.size()); i += 3)
		// 	{
		// 		const VEC3& p1 = position[ear_indices[i]];
		// 		const VEC3& p2 = position[ear_indices[i+1]];
		// 		const VEC3& p3 = position[ear_indices[i+2]];
		// 		if (intersection_ray_triangle(A, AB, p1, p2, p3, &intersection_point))
		// 		{
		// 			selected_per_thread[worker_index].push_back({ f, intersection_point, (intersection_point -
		// A).squaredNorm() }); 			i = uint32(ear_indices.size());
		// 		}
		// 	}
		// }
		return true;
	});

	std::vector<SelectedFace> result;

	for (const auto& selected_faces_vector : selected_per_thread)
		result.insert(result.end(), selected_faces_vector.begin(), selected_faces_vector.end());

	std::sort(result.begin(), result.end(),
			  [](const SelectedFace& f1, const SelectedFace& f2) -> bool { return std::get<2>(f1) < std::get<2>(f2); });

	return result;
}

} // namespace internal

template <typename MESH>
void picking(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position, const Vec3& A,
			 const Vec3& B, std::vector<typename mesh_traits<MESH>::Vertex>& result)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;
	using SelectedFace = std::tuple<Face, Vec3, Scalar>;

	std::vector<SelectedFace> selected_faces = internal::picking(m, vertex_position, A, B);

	CellMarkerStore<MESH, Vertex> cm(m);
	result.clear();
	result.reserve(selected_faces.size());
	for (const auto& sf : selected_faces)
	{
		Scalar min_d2 = std::numeric_limits<Scalar>::max();
		Vertex closest_vertex;

		Face f = std::get<0>(sf);
		const Vec3& I = std::get<1>(sf);

		foreach_incident_vertex(m, f, [&](Vertex v) -> bool {
			Scalar d2 = (value<Vec3>(m, vertex_position, v) - I).squaredNorm();
			if (d2 < min_d2)
			{
				min_d2 = d2;
				closest_vertex = v;
			}
			return true;
		});

		if (!cm.is_marked(closest_vertex))
		{
			cm.mark(closest_vertex);
			result.push_back(closest_vertex);
		}
	}
}

template <typename MESH>
void picking(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position, const Vec3& A,
			 const Vec3& B, std::vector<typename mesh_traits<MESH>::Edge>& result)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	using SelectedFace = std::tuple<Face, Vec3, Scalar>;

	std::vector<SelectedFace> selected_faces = internal::picking(m, vertex_position, A, B);

	CellMarkerStore<MESH, Edge> cm(m);
	result.clear();
	result.reserve(selected_faces.size());
	for (const auto& sf : selected_faces)
	{
		Scalar min_d2 = std::numeric_limits<Scalar>::max();
		Edge closest_edge;

		Face f = std::get<0>(sf);
		const Vec3& I = std::get<1>(sf);

		foreach_incident_edge(m, f, [&](Edge e) -> bool {
			std::vector<Vertex> vertices = incident_vertices(m, e);
			Scalar d2 = squared_distance_line_point(value<Vec3>(m, vertex_position, vertices[0]),
													value<Vec3>(m, vertex_position, vertices[1]), I);
			if (d2 < min_d2)
			{
				min_d2 = d2;
				closest_edge = e;
			}
			return true;
		});

		if (!cm.is_marked(closest_edge))
		{
			cm.mark(closest_edge);
			result.push_back(closest_edge);
		}
	}
}

template <typename MESH>
void picking(const MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position, const Vec3& A,
			 const Vec3& B, std::vector<typename mesh_traits<MESH>::Face>& result)
{
	using Face = typename mesh_traits<MESH>::Face;
	using SelectedFace = std::tuple<Face, Vec3, Scalar>;

	std::vector<SelectedFace> selected_faces = internal::picking(m, vertex_position, A, B);

	result.clear();
	result.reserve(selected_faces.size());
	for (const auto& sf : selected_faces)
		result.push_back(std::get<0>(sf));
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_PICKING_H_

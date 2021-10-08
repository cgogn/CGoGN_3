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

#ifndef CGOGN_GEOMETRY_ALGOS_MEDIAL_AXIS_H_
#define CGOGN_GEOMETRY_ALGOS_MEDIAL_AXIS_H_

#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <array>
#include <vector>

#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

namespace cgogn
{

namespace geometry
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;

inline Scalar compute_radius(const Vec3& p, const Vec3& n, const Vec3& q)
{
	// Compute radius of the ball that touches points p and q and whose center falls on the normal n from p
	Vec3 qp = p - q;
	Scalar d = qp.norm();
	// Scalar cos_theta = n.dot(p - q) / d;
	Scalar cos_theta = geometry::cos_angle(n, qp);
	return Scalar(d / (2 * cos_theta));
}

// adapted from https://github.com/tudelft3d/masbcpp

template <typename MESH>
void shrinking_ball_centers(MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
							const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal,
							typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_shrinking_ball_center)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	uint32 nb_vertices = nb_cells<Vertex>(m);
	uint32 nb_faces = nb_cells<Face>(m);

	auto bvh_vertex_index = add_attribute<uint32, Vertex>(m, "__bvh_vertex_index");
	std::vector<Vec3> bvh_vertex_position;
	bvh_vertex_position.reserve(nb_vertices);
	uint32 idx = 0;
	foreach_cell(m, [&](Vertex v) -> bool {
		value<uint32>(m, bvh_vertex_index, v) = idx++;
		bvh_vertex_position.push_back(value<Vec3>(m, vertex_position, v));
		return true;
	});
	std::vector<Face> bvh_faces;
	bvh_faces.reserve(nb_faces);
	std::vector<uint32> face_vertex_indices;
	face_vertex_indices.reserve(nb_faces * 3);
	foreach_cell(m, [&](Face f) -> bool {
		bvh_faces.push_back(f);
		foreach_incident_vertex(m, f, [&](Vertex v) -> bool {
			face_vertex_indices.push_back(value<uint32>(m, bvh_vertex_index, v));
			return true;
		});
		return true;
	});
	acc::BVHTree<uint32, Vec3>* surface_bvh = new acc::BVHTree<uint32, Vec3>(face_vertex_indices, bvh_vertex_position);

	std::vector<Vertex> kdt_vertices;
	std::vector<Vec3> kdt_vertex_position;
	kdt_vertices.reserve(nb_vertices);
	kdt_vertex_position.reserve(nb_vertices);
	foreach_cell(m, [&](Vertex v) -> bool {
		kdt_vertices.push_back(v);
		kdt_vertex_position.push_back(value<Vec3>(m, vertex_position, v));
		return true;
	});
	acc::KDTree<3, uint32>* surface_kdt = new acc::KDTree<3, uint32>(kdt_vertex_position);

	const Scalar denoise_preserve = 20.0 * M_PI / 180.0;
	const Scalar denoise_planar = 32.0 * M_PI / 180.0;
	const Scalar delta_convergence = 1e-5;
	const uint32 iteration_limit = 30;

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		const Vec3& p = value<Vec3>(m, vertex_position, v);
		const Vec3& n = value<Vec3>(m, vertex_normal, v);

		uint32 j = 0;
		Scalar r = 0.;

		acc::Ray<Vec3> ray{p, -n, 1e-10, acc::inf};
		acc::BVHTree<uint32, Vec3>::Hit h;
		if (surface_bvh->intersect(ray, &h))
		{
			Face f = bvh_faces[h.idx];
			std::vector<Vertex> vertices = incident_vertices(m, f);
			Vec3 ip = h.bcoords[0] * value<Vec3>(m, vertex_position, vertices[0]) +
					  h.bcoords[1] * value<Vec3>(m, vertex_position, vertices[1]) +
					  h.bcoords[2] * value<Vec3>(m, vertex_position, vertices[2]);
			r = (p - ip).norm() * 0.51;
		}
		// else
		// 	std::cout << "intersection point not found !!!";

		Vec3 c = p - (r * n);

		while (true)
		{
			// find closest point to c
			std::pair<uint32, Scalar> k_res;
			bool found = surface_kdt->find_nn(c, &k_res);
			// std::pair<uint32, Vec3> cp_res;
			// bool found = surface_bvh->closest_point(c, &cp_res);
			if (!found)
				std::cout << "closest point not found !!!";

			const Vec3& q = surface_kdt->vertex(k_res.first);
			Scalar d = k_res.second;
			// Vec3 q = cp_res.second;
			// Scalar d = (q - c).norm();

			// This should handle all (special) cases where we want to break the loop
			// - normal case when ball no longer shrinks
			// - the case where q == p
			// - any duplicate point cases
			if ((d >= r - delta_convergence) || (p == q))
				break;

			// Compute next ball center
			r = compute_radius(p, n, q);
			Vec3 c_next = p - (r * n);

			// // Denoising
			if (denoise_preserve > 0 || denoise_planar > 0)
			{
				Scalar separation_angle = geometry::angle(p - c_next, q - c_next);

				// if (j == 0 && denoise_planar > 0 && separation_angle < denoise_planar)
				// 	break;
				if (j > 0 && denoise_preserve > 0 && (separation_angle < denoise_preserve && r > (q - p).norm()))
					break;
			}

			// // Stop iteration if this looks like an infinite loop:
			if (j > iteration_limit)
				break;

			c = c_next;
			j++;
		}

		value<Vec3>(m, vertex_shrinking_ball_center, v) = c;

		return true;
	});

	delete surface_kdt;

	remove_attribute<Vertex>(m, bvh_vertex_index);
	delete surface_bvh;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_MEDIAL_AXIS_H_

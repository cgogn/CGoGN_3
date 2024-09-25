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

#include <cgogn/core/utils/numerics.h>
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

using geometry::Scalar;
using geometry::Vec3;

using namespace cgogn::numerics;

// Compute radius of the ball that touches points p and q and whose center falls on the normal n from p
inline Scalar compute_radius(const Vec3& p, const Vec3& n, const Vec3& q)
{
	Vec3 qp = p - q;
	Scalar d = qp.norm();
	// Scalar cos_theta = n.dot(p - q) / d;
	Scalar cos_theta = geometry::cos_angle(n, qp);
	return Scalar(d / (2 * cos_theta));
}

// const Scalar denoise_planar = 32.0 * M_PI / 180.0;
const Scalar denoise_preserve = 20.0 * M_PI / 180.0;
const Scalar delta_convergence = 1e-5;
const uint32 iteration_limit = 30;

template <typename MESH>
std::tuple<Vec3, Scalar, typename mesh_traits<MESH>::Vertex> shrinking_ball_center(
	const MESH& m, const Vec3& p, const Vec3& n,
	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
	const acc::BVHTree<uint32, Vec3>* surface_bvh, const std::vector<typename mesh_traits<MESH>::Face>& bvh_faces,
	const acc::KDTree<3, uint32>* surface_kdt, const std::vector<typename mesh_traits<MESH>::Vertex>& kdt_vertices,
	bool use_kdt_only = false, double initial_radius = 0.5)
{
	// initial radius is only used when use_kdt_only is true

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	uint32 j = 0;
	Scalar r = 0.0;

	// if (use_kdt_only)
		r = initial_radius;
	// else
	// {
	// 	acc::Ray<Vec3> ray{p, -n, 1e-10, acc::inf};
	// 	acc::BVHTree<uint32, Vec3>::Hit h;
	// 	if (surface_bvh->intersect(ray, &h))
	// 	{
	// 		Face f = bvh_faces[h.idx];
	// 		std::vector<Vertex> vertices = incident_vertices(m, f);
	// 		Vec3 ip = h.bcoords[0] * value<Vec3>(m, vertex_position, vertices[0]) +
	// 				  h.bcoords[1] * value<Vec3>(m, vertex_position, vertices[1]) +
	// 				  h.bcoords[2] * value<Vec3>(m, vertex_position, vertices[2]);
	// 		r = (p - ip).norm() * 0.75;
	// 	}
	// 	else
	// 		std::cout << "intersection point not found !!!";
	// }

	Vec3 c = p - (r * n);
	Vec3 q = p - (2 * r * n);
	Vertex q_v;

	while (true)
	{
		// Find closest point to c

		Scalar d;
		Vec3 q_next;
		Vertex q_next_v;

		// if (use_kdt_only)
		// {
		std::pair<uint32, Scalar> k_res;
		surface_kdt->find_nn(c, &k_res);
		q_next = surface_kdt->vertex(k_res.first);
		d = k_res.second;
		q_next_v = kdt_vertices[k_res.first];
		// }
		// else
		// {
		// 	std::pair<uint32, Vec3> cp_res;
		// 	surface_bvh->closest_point(c, &cp_res);
		// 	q_next = cp_res.second;
		// 	d = (q_next - c).norm();
		// 	Face f = bvh_faces[cp_res.first];
		// 	std::vector<Vertex> vertices = incident_vertices(m, f);
		// 	Scalar d0 = (q_next - value<Vec3>(m, vertex_position, vertices[0])).squaredNorm();
		// 	Scalar d1 = (q_next - value<Vec3>(m, vertex_position, vertices[1])).squaredNorm();
		// 	Scalar d2 = (q_next - value<Vec3>(m, vertex_position, vertices[2])).squaredNorm();
		// 	if (d0 < d1 && d0 < d2)
		// 		q_next_v = vertices[0];
		// 	else if (d1 < d0 && d1 < d2)
		// 		q_next_v = vertices[1];
		// 	else
		// 		q_next_v = vertices[2];
		// }

		// If the closest point is (almost) the same as the previous one, or if the ball no longer shrinks, we stop
		if (fabs(d - r) <= delta_convergence || (q_next - q).norm() < delta_convergence)
			break;

		// Compute next ball center
		Scalar r_next = compute_radius(p, n, q_next);
		Vec3 c_next = p - (r_next * n);

		// Denoising
		Scalar separation_angle = geometry::angle(p - c_next, q_next - c_next);
		if (j > 0 && separation_angle < denoise_preserve) // && r_next > // (q_next - p).norm())
			break;

		c = c_next;
		r = r_next;
		q = q_next;
		q_v = q_next_v;

		j++;
		if (j > iteration_limit)
			break;
	}

	return {c, r, q_v};
}

// template <typename MESH>
// std::tuple<Vec3, Scalar, Vec3> shrinking_ball_center(
// 	MESH& m, typename mesh_traits<MESH>::Vertex v,
// 	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
// 	const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal,
// 	const acc::BVHTree<uint32, Vec3>* surface_bvh, const std::vector<typename mesh_traits<MESH>::Face>& bvh_faces,
// 	const acc::KDTree<3, uint32>* surface_kdt)
// {
// 	using Vertex = typename mesh_traits<MESH>::Vertex;
// 	using Face = typename mesh_traits<MESH>::Face;

// 	const Vec3& p = value<Vec3>(m, vertex_position, v);
// 	const Vec3& n = value<Vec3>(m, vertex_normal, v);

// 	uint32 j = 0;
// 	Scalar r = 0.;

// 	acc::Ray<Vec3> ray{p, -n, 1e-5, acc::inf};
// 	acc::BVHTree<uint32, Vec3>::Hit h;
// 	if (surface_bvh->intersect(ray, &h))
// 	{
// 		Face f = bvh_faces[h.idx];
// 		std::vector<Vertex> vertices = incident_vertices(m, f);
// 		Vec3 ip = h.bcoords[0] * value<Vec3>(m, vertex_position, vertices[0]) +
// 				  h.bcoords[1] * value<Vec3>(m, vertex_position, vertices[1]) +
// 				  h.bcoords[2] * value<Vec3>(m, vertex_position, vertices[2]);
// 		r = (p - ip).norm() * 0.75;
// 	}
// 	// else
// 	// 	std::cout << "intersection point not found !!!";

// 	Vec3 c = p - (r * n);
// 	Vec3 q = p - (2 * r * n);

// 	while (true)
// 	{
// 		// find closest point to c
// 		std::pair<uint32, Scalar> k_res;
// 		bool found = surface_kdt->find_nn(c, &k_res);
// 		// std::pair<uint32, Vec3> cp_res;
// 		// bool found = surface_bvh->closest_point(c, &cp_res);
// 		if (!found)
// 			std::cout << "closest point not found !!!";

// 		const Vec3& q_next = surface_kdt->vertex(k_res.first);
// 		Scalar d = k_res.second;
// 		// Vec3 q_next = cp_res.second;
// 		// Scalar d = (q_next - c).norm();

// 		// This should handle all (special) cases where we want to break the loop
// 		// - normal case when ball no longer shrinks
// 		// - the case where q == p
// 		// - any duplicate point cases
// 		if ((d >= r - delta_convergence) || (p == q_next))
// 			break;

// 		// Compute next ball center
// 		r = compute_radius(p, n, q_next);
// 		Vec3 c_next = p - (r * n);

// 		// Denoising
// 		if (denoise_preserve > 0) // || denoise_planar > 0)
// 		{
// 			Scalar separation_angle = geometry::angle(p - c_next, q_next - c_next);

// 			// if (j == 0 && denoise_planar > 0 && separation_angle < denoise_planar)
// 			// 	break;
// 			if (j > 0 && denoise_preserve > 0 && (separation_angle < denoise_preserve && r > (q_next - p).norm()))
// 				break;
// 		}

// 		// Stop iteration if this looks like an infinite loop:
// 		if (j > iteration_limit)
// 			break;

// 		c = c_next;
// 		q = q_next;
// 		j++;
// 	}

// 	return {c, r, q};
// }

// adapted from https://github.com/tudelft3d/masbcpp

template <typename MESH>
void shrinking_ball_centers(MESH& m, const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position,
							const typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_normal,
							typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_shrinking_ball_center,
							typename mesh_traits<MESH>::template Attribute<Scalar>* vertex_shrinking_ball_radius,
							typename mesh_traits<MESH>::template Attribute<typename mesh_traits<MESH>::Vertex>*
								vertex_shrinking_ball_secondary_vertex)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	uint32 nb_vertices = nb_cells<Vertex>(m);
	uint32 nb_faces = nb_cells<Face>(m);

	auto bvh_vertex_index = add_attribute<uint32, Vertex>(m, "__bvh_vertex_index");

	std::vector<Vec3> vertex_position_vector;
	vertex_position_vector.reserve(nb_vertices);
	std::vector<Vertex> kdt_vertices;
	kdt_vertices.reserve(nb_vertices);
	uint32 idx = 0;
	foreach_cell(m, [&](Vertex v) -> bool {
		value<uint32>(m, bvh_vertex_index, v) = idx++;
		vertex_position_vector.push_back(value<Vec3>(m, vertex_position, v));
		kdt_vertices.push_back(v);
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

	acc::BVHTree<uint32, Vec3>* surface_bvh =
		new acc::BVHTree<uint32, Vec3>(face_vertex_indices, vertex_position_vector);

	acc::KDTree<3, uint32>* surface_kdt = new acc::KDTree<3, uint32>(vertex_position_vector);

	parallel_foreach_cell(m, [&](Vertex v) -> bool {
		auto [c, r, q] = shrinking_ball_center(m, value<Vec3>(m, vertex_position, v), value<Vec3>(m, vertex_normal, v),
											   vertex_position, surface_bvh, bvh_faces, surface_kdt, kdt_vertices);
		value<Vec3>(m, vertex_shrinking_ball_center, v) = c;
		value<Scalar>(m, vertex_shrinking_ball_radius, v) = r;
		value<Vertex>(m, vertex_shrinking_ball_secondary_vertex, v) = q;
		return true;
	});

	delete surface_kdt;

	remove_attribute<Vertex>(m, bvh_vertex_index);
	delete surface_bvh;
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_MEDIAL_AXIS_H_

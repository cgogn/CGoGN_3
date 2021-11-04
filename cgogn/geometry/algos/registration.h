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

#ifndef CGOGN_GEOMETRY_ALGOS_REGISTRATION_H_
#define CGOGN_GEOMETRY_ALGOS_REGISTRATION_H_

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/global.h>

#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/Sparse>
#include <libacc/bvh_tree.h>
#include <simpleICP/simpleicp.h>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
void rigid_register_mesh(MESH* source, typename mesh_traits<MESH>::template Attribute<Vec3>* source_vertex_position,
						 MESH* target,
						 const typename mesh_traits<MESH>::template Attribute<Vec3>* target_vertex_position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;

	uint32 source_nbv = nb_cells<Vertex>(*source);
	Eigen::MatrixXd X_source(source_nbv, 3);
	auto source_vertex_index = add_attribute<uint32, Vertex>(*source, "__vertex_index");
	uint32 source_vertex_idx = 0;
	foreach_cell(*source, [&](Vertex v) -> bool {
		const Vec3& p = value<Vec3>(*source, source_vertex_position, v);
		X_source(source_vertex_idx, 0) = p[0];
		X_source(source_vertex_idx, 1) = p[1];
		X_source(source_vertex_idx, 2) = p[2];
		value<uint32>(*source, source_vertex_index, v) = source_vertex_idx++;
		return true;
	});

	uint32 target_nbv = nb_cells<Vertex>(*target);
	Eigen::MatrixXd X_target(target_nbv, 3);
	auto target_vertex_index = add_attribute<uint32, Vertex>(*target, "__vertex_index");
	uint32 target_vertex_idx = 0;
	foreach_cell(*target, [&](Vertex v) -> bool {
		const Vec3& p = value<Vec3>(*target, target_vertex_position, v);
		X_target(target_vertex_idx, 0) = p[0];
		X_target(target_vertex_idx, 1) = p[1];
		X_target(target_vertex_idx, 2) = p[2];
		value<uint32>(*target, target_vertex_index, v) = target_vertex_idx++;
		return true;
	});

	Mat4 t = SimpleICP(X_target, X_source);

	Eigen::MatrixXd X_sourceH(4, source_nbv);
	foreach_cell(*source, [&](Vertex v) -> bool {
		uint32 idx = value<uint32>(*source, source_vertex_index, v);
		const Vec3& p = value<Vec3>(*source, source_vertex_position, v);
		X_sourceH(0, idx) = p[0];
		X_sourceH(1, idx) = p[1];
		X_sourceH(2, idx) = p[2];
		X_sourceH(3, idx) = 1.0;
		return true;
	});
	Eigen::MatrixXd res = t * X_sourceH;
	foreach_cell(*source, [&](Vertex v) -> bool {
		uint32 idx = value<uint32>(*source, source_vertex_index, v);
		Vec3& p = value<Vec3>(*source, source_vertex_position, v);
		p[0] = res(0, idx);
		p[1] = res(1, idx);
		p[2] = res(2, idx);
		return true;
	});

	remove_attribute<Vertex>(*source, source_vertex_index);
	remove_attribute<Vertex>(*target, target_vertex_index);
}

template <typename MESH>
void non_rigid_register_mesh(MESH* source, typename mesh_traits<MESH>::template Attribute<Vec3>* source_vertex_position,
							 MESH* target,
							 const typename mesh_traits<MESH>::template Attribute<Vec3>* target_vertex_position)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	// index source vertices
	auto source_vertex_index = add_attribute<uint32, Vertex>(*source, "vertex_index");
	uint32 source_nb_vertices = 0;
	foreach_cell(*source, [&](Vertex v) -> bool {
		value<uint32>(*source, source_vertex_index, v) = source_nb_vertices++;
		return true;
	});

	// index target vertices
	// & build BVH
	auto target_vertex_index = add_attribute<uint32, Vertex>(*target, "vertex_index");
	std::vector<Vec3> target_vertex_positions;
	uint32 target_nb_vertices = 0;
	foreach_cell(*target, [&](Vertex v) -> bool {
		value<uint32>(*target, target_vertex_index, v) = target_nb_vertices++;
		target_vertex_positions.push_back(value<Vec3>(*target, target_vertex_position, v));
		return true;
	});
	// std::vector<Face> target_faces;
	std::vector<uint32> target_face_vertex_indices;
	foreach_cell(*target, [&](Face f) -> bool {
		// target_faces.push_back(f);
		foreach_incident_vertex(*target, f, [&](Vertex v) -> bool {
			target_face_vertex_indices.push_back(value<uint32>(*target, target_vertex_index, v));
			return true;
		});
		return true;
	});
	acc::BVHTree<uint32, Vec3> target_bvh(target_face_vertex_indices, target_vertex_positions);

	// compute source topo laplacian
	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> LAPL =
		geometry::topo_laplacian_matrix(*source, source_vertex_index.get());
	// LAPL = LAPL.diagonal().asDiagonal().inverse() * LAPL;

	Eigen::MatrixXd vpos(source_nb_vertices, 3);
	Eigen::MatrixXd lapl(source_nb_vertices, 3);

	parallel_foreach_cell(*source, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(*source, source_vertex_index, v);
		const Vec3& pos = value<Vec3>(*source, source_vertex_position, v);
		vpos(vidx, 0) = pos[0];
		vpos(vidx, 1) = pos[1];
		vpos(vidx, 2) = pos[2];
		return true;
	});

	lapl = LAPL * vpos;

	// setup matrix
	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(2 * source_nb_vertices, source_nb_vertices);
	std::vector<Eigen::Triplet<Scalar>> Acoeffs;
	Acoeffs.reserve(source_nb_vertices * 10);
	Eigen::MatrixXd b(2 * source_nb_vertices, 3);

	Scalar fit = 0.1;

	foreach_cell(*source, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(*source, source_vertex_index, v);
		uint32 nbv = 0;
		foreach_adjacent_vertex_through_edge(*source, v, [&](Vertex av) -> bool {
			uint32 avidx = value<uint32>(*source, source_vertex_index, av);
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(avidx), 1));
			++nbv;
			return true;
		});
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(vidx), -1 * Scalar(nbv)));
		Acoeffs.push_back(Eigen::Triplet<Scalar>(int(source_nb_vertices + vidx), int(vidx), fit));
		b(vidx, 0) = lapl(vidx, 0);
		b(vidx, 1) = lapl(vidx, 1);
		b(vidx, 2) = lapl(vidx, 2);
		Vec3 pos = target_bvh.closest_point(value<Vec3>(*source, source_vertex_position, v));
		b(source_nb_vertices + vidx, 0) = fit * pos[0];
		b(source_nb_vertices + vidx, 1) = fit * pos[1];
		b(source_nb_vertices + vidx, 2) = fit * pos[2];
		return true;
	});
	A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

	Eigen::SparseMatrix<Scalar, Eigen::ColMajor> At = A.transpose();
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(At * A);

	vpos = solver.solve(At * b);

	parallel_foreach_cell(*source, [&](Vertex v) -> bool {
		uint32 vidx = value<uint32>(*source, source_vertex_index, v);
		Vec3& pos = value<Vec3>(*source, source_vertex_position, v);
		pos[0] = vpos(vidx, 0);
		pos[1] = vpos(vidx, 1);
		pos[2] = vpos(vidx, 2);
		return true;
	});

	remove_attribute<Vertex>(*source, source_vertex_index);
	remove_attribute<Vertex>(*target, target_vertex_index);
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_REGISTRATION_H_

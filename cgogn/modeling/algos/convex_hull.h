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

#ifndef CGOGN_MODELING_ALGOS_CONVEX_HULL_H_
#define CGOGN_MODELING_ALGOS_CONVEX_HULL_H_

#include <cgogn/core/utils/numerics.h>

#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_ops/volume.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/functions/distance.h>
#include <cgogn/geometry/functions/normal.h>
#include <cgogn/geometry/functions/orientation.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <array>
#include <vector>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;

///////////
// CMap2 //
///////////

std::array<CMap2::Vertex, 4> tet_vertices(CMap2& m, CMap2::Volume v)
{
	return {CMap2::Vertex(v.dart), CMap2::Vertex(phi1(m, v.dart)), CMap2::Vertex(phi_1(m, v.dart)),
			CMap2::Vertex(phi<-1, 2, -1>(m, v.dart))};
}

// return [horizon_halfedges, visible_faces]
std::pair<std::vector<CMap2::HalfEdge>, std::vector<CMap2::Face>> build_horizon(CMap2& m,
																				CMap2::Attribute<Vec3>* vertex_position,
																				const Vec3& point, CMap2::Face f)
{
	std::vector<CMap2::HalfEdge> horizon_halfedges;
	std::vector<CMap2::Face> visible_faces;
	CellMarkerStore<CMap2, CMap2::Face> visible(m);
	visible_faces.push_back(f);
	visible.mark(f);
	for (uint32 i = 0; i < uint32(visible_faces.size()); ++i)
	{
		const CMap2::Face f_i = visible_faces[i];
		foreach_adjacent_face_through_edge(m, f_i, [&](CMap2::Face af) -> bool {
			if (visible.is_marked(af))
				return true;
			std::vector<CMap2::Vertex> vertices = incident_vertices(m, af);
			Vec3 N = geometry::normal(value<Vec3>(m, vertex_position, vertices[0]),
									  value<Vec3>(m, vertex_position, vertices[1]),
									  value<Vec3>(m, vertex_position, vertices[2]));
			Scalar d = -N.dot(value<Vec3>(m, vertex_position, vertices[0]));
			Scalar dist = geometry::signed_distance_plane_point(N, d, point);
			if (dist > 0)
			{
				visible_faces.push_back(af);
				visible.mark(af);
			}
			else
				horizon_halfedges.push_back(CMap2::HalfEdge(phi2(m, af.dart)));
			return true;
		});
	}

	// reorder horizon halfedges into a cycle
	for (uint32 i = 0; i < uint32(horizon_halfedges.size()) - 1; ++i)
	{
		uint32 end_vertex_index = index_of(m, CMap2::Vertex(phi1(m, horizon_halfedges[i].dart)));
		for (uint32 j = i + 1; j < horizon_halfedges.size(); ++j)
		{
			uint32 begin_vertex_index = index_of(m, CMap2::Vertex(horizon_halfedges[i].dart));
			if (begin_vertex_index == end_vertex_index)
			{
				if (j > i + 1)
					std::swap(horizon_halfedges[i + 1], horizon_halfedges[j]);
				break;
			}
		}
	}

	return {horizon_halfedges, visible_faces};
}

CMap2::Vertex remove_visible_faces_and_fill(CMap2& m, std::vector<CMap2::HalfEdge>& horizon_halfedges,
											std::vector<CMap2::Face>& visible_faces)
{
	Dart d = phi2(m, horizon_halfedges[0].dart);
	for (CMap2::HalfEdge he : horizon_halfedges)
		phi2_unsew(m, he.dart);
	for (CMap2::Face f : visible_faces)
		remove_face(static_cast<CMap1&>(m), CMap1::Face(f.dart));

	CMap2::Face f = close_hole(m, d);
	if (is_indexed<CMap2::Face>(m))
	{
		if (index_of(m, f) == INVALID_INDEX)
			set_index(m, f, new_index<CMap2::Face>(m));
	}

	Dart d0 = phi1(m, f.dart);

	cut_face(m, CMap2::Vertex(d0), CMap2::Vertex(f.dart));
	cut_edge(m, CMap2::Edge(phi_1(m, d0)));

	Dart x = phi<-1, -1>(m, d0);
	Dart dd = phi1(m, d0);
	while (dd != x)
	{
		cut_face(m, CMap2::Vertex(dd), CMap2::Vertex(phi1(m, x)));
		dd = phi1(m, dd);
	}

	return CMap2::Vertex(phi2(m, x));
}

/////////////
// GENERIC //
/////////////

// adapted from https://github.com/akuukka/quickhull

template <typename MESH>
void convex_hull(const std::vector<Vec3>& points, MESH& m,
				 typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using HalfEdge = typename mesh_traits<MESH>::HalfEdge;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;
	using Volume = typename mesh_traits<MESH>::Volume;

	if (points.size() < 4)
		return;

	// if only 4 points
	if (points.size() == 4)
	{
		std::array<uint32, 4> ids = {0, 1, 2, 3};
		if (geometry::test_orientation_3D(points[0], points[1], points[2], points[3]) == geometry::Orientation3D::UNDER)
			ids = {0, 2, 1, 3};

		Volume tet = add_pyramid(m, 3);
		std::array<Vertex, 4> vertices = tet_vertices(m, tet);

		value<Vec3>(m, vertex_position, vertices[0]) = points[ids[0]];
		value<Vec3>(m, vertex_position, vertices[1]) = points[ids[1]];
		value<Vec3>(m, vertex_position, vertices[2]) = points[ids[2]];
		value<Vec3>(m, vertex_position, vertices[3]) = points[ids[3]];

		return;
	}

	// compute points bounding box
	// and keep index of extreme points (xmax, xmin, ymax, ymin, zmax, zmin)
	Vec3 bb_min, bb_max;
	std::array<uint32, 6> extreme_points_index = {0, 0, 0, 0, 0, 0};
	for (uint32 i = 0; i < 3; ++i)
	{
		bb_min[i] = std::numeric_limits<Scalar>::max();
		bb_max[i] = std::numeric_limits<Scalar>::lowest();
	}
	for (uint32 i = 0; i < points.size(); ++i)
	{
		if (points[i][0] < bb_min[0])
		{
			bb_min[0] = points[i][0];
			extreme_points_index[0] = i;
		}
		if (points[i][0] > bb_max[0])
		{
			bb_max[0] = points[i][0];
			extreme_points_index[1] = i;
		}
		if (points[i][1] < bb_min[1])
		{
			bb_min[1] = points[i][1];
			extreme_points_index[2] = i;
		}
		if (points[i][1] > bb_max[1])
		{
			bb_max[1] = points[i][1];
			extreme_points_index[3] = i;
		}
		if (points[i][2] < bb_min[2])
		{
			bb_min[2] = points[i][2];
			extreme_points_index[4] = i;
		}
		if (points[i][2] > bb_max[2])
		{
			bb_max[2] = points[i][2];
			extreme_points_index[5] = i;
		}
	}

	// compute scale (bigger bb component value)
	Scalar scale = 0;
	for (uint32 i = 0; i < 3; ++i)
	{
		Scalar av = std::max(std::abs(bb_min[i]), std::abs(bb_max[i]));
		if (av > scale)
			scale = av;
	}

	Scalar epsilon = 0.0000001 * scale;
	Scalar epsilon_squared = epsilon * epsilon;

	// find 2 most distant extreme points
	Scalar max_d = epsilon_squared;
	std::pair<uint32, uint32> selected_points_index;
	for (uint32 i = 0; i < 6; ++i)
	{
		for (uint32 j = i + 1; j < 6; ++j)
		{
			Scalar d = (points[extreme_points_index[i]] - points[extreme_points_index[j]]).squaredNorm();
			if (d > max_d)
			{
				max_d = d;
				selected_points_index = {extreme_points_index[i], extreme_points_index[j]};
			}
		}
	}
	if (max_d == epsilon_squared)
		return;
	assert(selected_points_index.first != selected_points_index.second);

	// find the most distant point to the line between the 2 chosen extreme points
	max_d = epsilon_squared;
	uint32 most_distant_index = 0;
	for (uint32 i = 0; i < points.size(); ++i)
	{
		Scalar dist = geometry::squared_distance_line_point(points[selected_points_index.first],
															points[selected_points_index.second], points[i]);
		if (dist > max_d)
		{
			max_d = dist;
			most_distant_index = i;
		}
	}
	if (max_d == epsilon_squared)
		return;
	assert(selected_points_index.first != most_distant_index && selected_points_index.second != most_distant_index);

	// these three points form the base triangle for our tet
	std::array<uint32, 3> base_triangle_index{selected_points_index.first, selected_points_index.second,
											  most_distant_index};
	Vec3 base_triangle_pos[] = {points[base_triangle_index[0]], points[base_triangle_index[1]],
								points[base_triangle_index[2]]};

	// next step is to find the 4th vertex of the tetrahedron
	// we naturally choose the point farthest away from the triangle plane
	max_d = epsilon;
	most_distant_index = 0;
	Vec3 N = geometry::normal(base_triangle_pos[0], base_triangle_pos[1], base_triangle_pos[2]);
	Scalar d = -(base_triangle_pos[0].dot(N));
	for (uint32 i = 0; i < points.size(); ++i)
	{
		Scalar dist = geometry::distance_plane_point(N, d, points[i]);
		if (dist > max_d)
		{
			max_d = dist;
			most_distant_index = i;
		}
	}
	if (max_d == epsilon_squared)
		return;

	// Create the initial tetrahedron from the selected points
	std::array<uint32, 4> ids = {0, 1, 2, 3};
	if (geometry::test_orientation_3D(base_triangle_pos[0], base_triangle_pos[1], base_triangle_pos[2],
									  points[most_distant_index]) == geometry::Orientation3D::UNDER)
		ids = {0, 2, 1, 3};

	Volume tet = add_pyramid(m, 3);
	std::array<Vertex, 4> vertices = tet_vertices(m, tet);

	value<Vec3>(m, vertex_position, vertices[0]) = base_triangle_pos[ids[0]];
	value<Vec3>(m, vertex_position, vertices[1]) = base_triangle_pos[ids[1]];
	value<Vec3>(m, vertex_position, vertices[2]) = base_triangle_pos[ids[2]];
	value<Vec3>(m, vertex_position, vertices[3]) = points[most_distant_index];

	auto face_points_on_positive_side = add_attribute<std::vector<uint32>, Face>(m, "points_on_positive_side");
	auto face_most_distant_point_dist = add_attribute<Scalar, Face>(m, "most_distant_point_dist");
	face_most_distant_point_dist->fill(0);
	auto face_most_distant_point_index = add_attribute<uint32, Face>(m, "most_distant_point_index");
	face_most_distant_point_index->fill(0);

	// register points outside the tetrahedron in the faces
	for (uint32 i = 0; i < points.size(); ++i)
	{
		foreach_cell(m, [&](Face f) -> bool {
			std::vector<Vertex> vertices = incident_vertices(m, f);
			Vec3 N = geometry::normal(value<Vec3>(m, vertex_position, vertices[0]),
									  value<Vec3>(m, vertex_position, vertices[1]),
									  value<Vec3>(m, vertex_position, vertices[2]));
			Scalar d = -N.dot(value<Vec3>(m, vertex_position, vertices[0]));
			Scalar dist = geometry::signed_distance_plane_point(N, d, points[i]);
			if (dist > 0) // && dist * dist > epsilon_squared * N.squaredNorm())
			{
				value<std::vector<uint32>>(m, face_points_on_positive_side, f).push_back(i);
				if (dist > value<Scalar>(m, face_most_distant_point_dist, f))
				{
					value<Scalar>(m, face_most_distant_point_dist, f) = dist;
					value<uint32>(m, face_most_distant_point_index, f) = i;
				}
				return false;
			}
			return true;
		});
	}

	// init face stack with faces that have points assigned to them
	std::vector<Face> face_list;
	foreach_cell(m, [&](Face f) -> bool {
		if (value<std::vector<uint32>>(m, face_points_on_positive_side, f).size() > 0)
			face_list.push_back(f);
		return true;
	});

	while (!face_list.empty())
	{
		Face f = face_list.back();
		face_list.pop_back();

		// pick the most distant point to this triangle plane as the point to which we extrude
		const uint32 active_point_index = value<uint32>(m, face_most_distant_point_index, f);
		const Vec3& active_point = points[active_point_index];

		// create the list of horizon halfedges
		auto [horizon_halfedges, visible_faces] = build_horizon(m, vertex_position, active_point, f);

		// save visible faces points
		std::vector<uint32> visible_points;
		for (Face f : visible_faces)
		{
			const auto& vp = value<std::vector<uint32>>(m, face_points_on_positive_side, f);
			visible_points.insert(visible_points.end(), vp.begin(), vp.end());
		}

		// remove faces & fill hole with new triangles
		Vertex v = remove_visible_faces_and_fill(m, horizon_halfedges, visible_faces);
		value<Vec3>(m, vertex_position, v) = active_point;

		foreach_incident_face(m, v, [&](Face iface) -> bool {
			value<std::vector<uint32>>(m, face_points_on_positive_side, iface).clear();
			value<Scalar>(m, face_most_distant_point_dist, iface) = 0;
			value<uint32>(m, face_most_distant_point_index, iface) = 0;
			return true;
		});

		// register points in the new faces
		for (uint32 i : visible_points)
		{
			if (i == active_point_index)
				continue;
			foreach_incident_face(m, v, [&](Face iface) -> bool {
				std::vector<Vertex> vertices = incident_vertices(m, iface);
				Vec3 N = geometry::normal(value<Vec3>(m, vertex_position, vertices[0]),
										  value<Vec3>(m, vertex_position, vertices[1]),
										  value<Vec3>(m, vertex_position, vertices[2]));
				Scalar d = -N.dot(value<Vec3>(m, vertex_position, vertices[0]));
				Scalar dist = geometry::signed_distance_plane_point(N, d, points[i]);
				if (dist > 0) // && dist * dist > epsilon_squared * N.squaredNorm())
				{
					value<std::vector<uint32>>(m, face_points_on_positive_side, iface).push_back(i);
					if (dist > value<Scalar>(m, face_most_distant_point_dist, iface))
					{
						value<Scalar>(m, face_most_distant_point_dist, iface) = dist;
						value<uint32>(m, face_most_distant_point_index, iface) = i;
					}
					return false;
				}
				return true;
			});
		}

		// add faces in stack
		foreach_incident_face(m, v, [&](Face iface) -> bool {
			if (value<std::vector<uint32>>(m, face_points_on_positive_side, iface).size() > 0)
				face_list.push_back(iface);
			return true;
		});
	}

	remove_attribute<Face>(m, face_points_on_positive_side);
	remove_attribute<Face>(m, face_most_distant_point_dist);
	remove_attribute<Face>(m, face_most_distant_point_index);
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_CONVEX_HULL_H_

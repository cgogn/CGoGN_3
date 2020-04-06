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

#ifndef CGOGN_GEOMETRY_ALGOS_GRAPH_TO_HEX_H_
#define CGOGN_GEOMETRY_ALGOS_GRAPH_TO_HEX_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/io/surface/surface_import.h>
namespace cgogn
{

template <typename MESH, typename CELL>
class CellMarker;

namespace modeling
{

using Vec3 = geometry::Vec3;
using Vec3i = geometry::Vec3i;
using Scalar = geometry::Scalar;
using Mat3 = geometry::Mat3;

struct GData
{
	std::vector<std::pair<Graph::HalfEdge, Graph::HalfEdge>> branches;
	std::vector<Graph::Vertex> intersections;
};

struct GAttributes
{
	std::shared_ptr<Graph::Attribute<Vec3>> vertex_position;
	std::shared_ptr<Graph::Attribute<Scalar>> vertex_radius;
	std::shared_ptr<Graph::Attribute<Dart>> vertex_contact_surface;
	std::shared_ptr<Graph::Attribute<Dart>> halfedge_contact_surface_face;
	std::shared_ptr<Graph::Attribute<Mat3>> halfedge_frame;
};

struct M2Attributes
{
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_position;
	std::shared_ptr<CMap2::Attribute<Dart>> dual_vertex_graph_branch;
	std::shared_ptr<CMap2::Attribute<Graph::Vertex>> volume_gvertex;
	std::shared_ptr<CMap2::Attribute<Vec3>> volume_center;
	std::shared_ptr<CMap2::Attribute<Vec3>> edge_mid;
	std::shared_ptr<CMap2::Attribute<Dart>> halfedge_volume_connection;
};

struct M3Attributes
{
	std::shared_ptr<CMap3::Attribute<Vec3>> vertex_position;
	std::shared_ptr<CMap3::Attribute<Mat3>> corner_frame;
	std::shared_ptr<CMap3::Attribute<Mat3>> hex_frame;

	std::shared_ptr<CMap3::Attribute<Scalar>> scaled_jacobian;
	std::shared_ptr<CMap3::Attribute<Scalar>> jacobian;
	std::shared_ptr<CMap3::Attribute<Scalar>> max_frobenius;
	std::shared_ptr<CMap3::Attribute<Scalar>> mean_frobenius;

	std::shared_ptr<CMap3::Attribute<Vec3>> color_scaled_jacobian;
	std::shared_ptr<CMap3::Attribute<Vec3>> color_jacobian;
	std::shared_ptr<CMap3::Attribute<Vec3>> color_max_frobenius;
	std::shared_ptr<CMap3::Attribute<Vec3>> color_mean_frobenius;
};

bool graph_to_hex(Graph& g, CMap2& m2, CMap3& m3);

/*****************************************************************************/
/* utils                                                                     */
/*****************************************************************************/

void index_volume_cells(CMap2& m, CMap2::Volume vol);
void sew_volumes(CMap3& m, Dart d0, Dart d1);
Graph::HalfEdge branch_extremity(const Graph& g, Graph::HalfEdge h, CellMarker<Graph, Graph::Edge>& cm);
Dart add_branch_section(CMap3& m3);
void project_on_sphere(Vec3& P, const Vec3& C, Scalar R);
void shift_frame(Mat3& frame, uint32 nb_shifts);
void dualize_volume(CMap2& m, CMap2::Volume vol, M2Attributes& m2Attribs, const Graph& g, GAttributes& gAttribs);
bool dijkstra_topo(CMap2& m2, CMap2::Vertex v0, std::shared_ptr<CMap2::Attribute<Dart>> previous,
				   std::shared_ptr<CMap2::Attribute<uint32>> dist);
Dart convex_hull(CMap2& m2, const cgogn::io::SurfaceImportData& surface_data);
Dart remesh(CMap2& m, CMap2::Volume vol, M2Attributes& m2Attribs);
Vec3 slerp(Vec3 A, Vec3 B, Scalar coef, bool in);
Scalar angle_on_sphere(Vec3 A, Vec3 B, Vec3 C);
Scalar edge_max_angle(CMap2& m2, CMap2::Edge e, M2Attributes& m2Attribs);
Scalar min_cut_angle(CMap2& m2, CMap2::Vertex v0, CMap2::Vertex v1, M2Attributes& m2Attribs);
Vec3 spherical_barycenter(std::vector<Vec3>& points, uint32 iterations);
void extract_volume_surface(CMap3& m3, CMap2& m2);
void catmull_clark_approx(CMap2& m2, uint32 iterations);
void catmull_clark_inter(CMap2& m, uint32 iterations);


/*****************************************************************************/
/* data preparation                                                          */
/*****************************************************************************/
bool subdivide_graph(Graph& g);
bool get_graph_data(const Graph& g, GData& gData);
bool add_graph_attributes(Graph& g, GAttributes& gAttribs);
bool add_cmap2_attributes(CMap2& m2, M2Attributes& m2Attribs);

/*****************************************************************************/
/* contact surfaces generation                                               */
/*****************************************************************************/

bool build_contact_surfaces(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs);
void build_contact_surface_1(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
							 Graph::Vertex v);
void build_contact_surface_2(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
							 Graph::Vertex v);
void build_contact_surface_3(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
							 Graph::Vertex v);
void build_contact_surface_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
							 Graph::Vertex v);

/*****************************************************************************/
/* frames initialization & propagation                                       */
/*****************************************************************************/

bool create_intersection_frames(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs);
bool create_intersection_frame_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
								 Graph::Vertex v);

bool propagate_frames(const Graph& g, GAttributes& gAttribs, const GData& gData, CMap2& m2);
void propagate_frame_n_1(const Graph& g, GAttributes& gAttribs, Graph::HalfEdge h_from_start);
bool propagate_frame_n_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, Graph::HalfEdge h_from_start);

/*****************************************************************************/
/* contact surfaces geometry                                                 */
/*****************************************************************************/

bool set_contact_surfaces_geometry(const Graph& g, const GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs);

/*****************************************************************************/
/* volume mesh generation                                                    */
/*****************************************************************************/

bool build_branch_sections(Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3);
bool sew_branch_sections(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3);
bool set_volumes_geometry(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3);

/*****************************************************************************/
/* mesh volume quality                                                       */
/*****************************************************************************/

bool add_quality_attributes(CMap3& m3, M3Attributes& m3Attribs);
bool set_hex_frames(CMap3& m3, M3Attributes& m3Attribs);
bool compute_scaled_jacobians(CMap3& m3, M3Attributes& m3Attribs, bool add_color);
bool compute_jacobians(CMap3& m3, M3Attributes& m3Attribs, bool add_color);
bool compute_maximum_aspect_frobenius(CMap3& m3, M3Attributes& m3Attribs, bool add_color);
bool compute_mean_aspect_frobenius(CMap3& m3, M3Attributes& m3Attribs, bool add_color);
Scalar frame_frobenius(Mat3 frame);
Vec3 get_quality_color(Scalar quality);



} // namespace modeling
} // namespace cgogn
#endif // CGOGN_GEOMETRY_ALGOS_GRAPH_TO_HEX_H_

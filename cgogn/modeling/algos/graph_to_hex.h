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

#ifndef CGOGN_MODELING_ALGOS_GRAPH_TO_HEX_H_
#define CGOGN_MODELING_ALGOS_GRAPH_TO_HEX_H_

#include <cgogn/core/types/cells_set.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/io/surface/surface_import.h>
#include <cgogn/modeling/algos/graph_utils.h>

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

struct GAttributes
{
	std::shared_ptr<Graph::Attribute<Vec3>> vertex_position;
	std::shared_ptr<Graph::Attribute<Scalar>> vertex_radius;
	std::shared_ptr<Graph::Attribute<Dart>> vertex_contact_surface;
	std::shared_ptr<Graph::Attribute<Dart>> halfedge_volume_connection;
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
	std::shared_ptr<CMap2::Attribute<CMap2*>> ortho_scaffold;
};

struct M3Attributes
{
	std::shared_ptr<CMap3::Attribute<Vec3>> vertex_position;
	std::shared_ptr<CMap3::Attribute<Graph::HalfEdge>> volume_graph_connection;

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

std::tuple<GAttributes, M2Attributes, M3Attributes> graph_to_hex(Graph& g, CMap2& m2, CMap3& m3);

/*****************************************************************************/
/* utils                                                                     */
/*****************************************************************************/

void index_volume_cells(CMap2& m, CMap2::Volume vol);
void sew_volumes(CMap3& m, Dart d0, Dart d1);
void unsew_volumes(CMap3& m, Dart d0);
Dart add_branch_section(CMap3& m3);
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

void catmull_clark_approx(CMap2& m2, CMap2::Attribute<Vec3>* vertex_position, uint32 iterations);
void catmull_clark_inter(CMap2& m, CMap2::Attribute<Vec3>* vertex_position, uint32 iterations);

// void bloat(CMap3& m3);
void bloat(CMap3& m3, const Graph& g, const GAttributes& gAttribs);

void padding(CMap3& m3);

void create_ortho_hex(const Graph& g, CMap2& m2, CMap2& contact_surface, CMap3& m3);
bool find_frame(const Graph& g, Graph::Vertex gv, Mat3& frame);
bool find_inter_frame(const Graph& g, Graph::Vertex gv, const GAttributes& gAttribs, Mat3& frame);

/*****************************************************************************/
/* subdivision                                                               */
/*****************************************************************************/

// void mark_tranversal_faces(CMap3& m3, CMap2& m2, M2Attributes& m2Attribs, CellMarker<CMap3, CMap3::Face>& cm);

// void subdivide_length_wise(CMap3& m3, M3Attributes& m3Attribs, CellMarker<CMap3, CMap3::Face>& trans_faces, Graph& g,
// 						   GAttributes& gAttribs);
// void subdivide_width_wise(CMap3& m3, M3Attributes& m3Attribs, CellMarker<CMap3, CMap3::Face>& trans_faces, Graph& g,
// 						  GAttributes& gAttribs);
// void trisect_length_wise(CMap3& m3, M3Attributes& m3Attribs, CellMarker<CMap3, CMap3::Face>& trans_faces, Graph& g,
// 						 GAttributes& gAttribs);
// void get_loop_path(CMap3& m3, Dart d0, std::vector<Dart>& path);
// void cut_chunk(CMap3& m3, M3Attributes& m3Attribs, CellMarker<CMap3, CMap3::Face>& trans_faces, Graph& g,
// 			   GAttributes& gAttribs, Graph::Edge eg, Scalar slice);

CMap3::Edge find_fiber_dir(CMap3& m3, CMap3::Face f);
uint32 get_ring_size(CMap3& m3, CMap3::Edge e);
bool unchecked_ring(CMap3& m3, CMap3::Edge e, uint32 ring_size, CellMarker<CMap3, CMap3::Edge>& visited_edge);
void cut_slice(CMap3& m3, CMap3::Attribute<Vec3>* vertex_position, CellCache<CMap3>& slice);
CellCache<CMap3> get_slice(CMap3& m, CMap3::Edge e);
void volume_fiber_spread(CMap3& m, CellCache<CMap3>& surface_fibers, CellMarker<CMap3, CMap3::Edge>& edge_fibers);
CellCache<CMap3> surface_fiber_spread(CMap3& m, CMap3::Edge e0);
void mark_mesh_fibers(CMap3& m3, CMap3::Edge e, CellMarker<CMap3, CMap3::Edge>& edge_fibers);
void fiber_aligned_subdivision(CMap3& m, CellMarker<CMap3, CMap3::Edge>& fibers);


// void mark_fibers_from_input();
// void fiber_aligned_subdivision(CMap3& m3);

/*****************************************************************************/
/* data preparation                                                          */
/*****************************************************************************/

bool subdivide_graph(Graph& g);
bool add_graph_attributes(Graph& g, GAttributes& gAttribs);
bool add_cmap2_attributes(CMap2& m2, M2Attributes& m2Attribs);
bool add_cmap3_attributes(CMap3& m3, M3Attributes& m3Attribs);

/*****************************************************************************/
/* contact surfaces generation                                               */
/*****************************************************************************/

bool build_contact_surfaces(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs);
void build_contact_surface_1(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
							 Graph::Vertex v);
void build_contact_surface_2(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
							 Graph::Vertex v);
void build_contact_surface_orange(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
								  Graph::Vertex v);
// void build_contact_surface_3(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
// 							 Graph::Vertex v);
void build_contact_surface_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
							 Graph::Vertex v);
bool build_contact_surface_ortho(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
								 Graph::Vertex v);

/*****************************************************************************/
/* frames initialization & propagation                                       */
/*****************************************************************************/

bool create_intersection_frames(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs);
bool create_intersection_frame_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
								 Graph::Vertex v);
bool create_extremity_frame(const Graph& g, GAttributes& gAttribs, Graph::Vertex v);

bool propagate_frames(const Graph& g, GAttributes& gAttribs, const GraphData& gData, CMap2& m2);
void propagate_frame_n_1(const Graph& g, GAttributes& gAttribs, Graph::HalfEdge h_from_start);
bool propagate_frame_n_n(const Graph& g, GAttributes& gAttribs, CMap2& m2, Graph::HalfEdge h_from_start);

/*****************************************************************************/
/* contact surfaces geometry                                                 */
/*****************************************************************************/

bool set_contact_surfaces_geometry(const Graph& g, const GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs);
bool set_contact_surfaces_geometry_from_surface(const Graph& g, const GAttributes& gAttribs, CMap2& m2,
												M2Attributes& m2Attribs, const CMap2& surface);

/*****************************************************************************/
/* volume mesh generation                                                    */
/*****************************************************************************/
void insert_ortho_chunks(Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3);
bool build_branch_sections(Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3);
bool sew_branch_sections(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3);
bool set_volumes_geometry(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3, M3Attributes& m3Attribs);
// bool set_volumes_geometry(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3);

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_GRAPH_TO_HEX_H_

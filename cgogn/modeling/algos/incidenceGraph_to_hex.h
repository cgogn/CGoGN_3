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

#ifndef CGOGN_MODELING_ALGOS_INCIDENCEGRAPH_TO_HEX_H_
#define CGOGN_MODELING_ALGOS_INCIDENCEGRAPH_TO_HEX_H_

#include <cgogn/core/types/cmap/cmap3.h>
#include <cgogn/core/types/incidence_graph/incidence_graph.h>
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/graph_utils.h>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;
using Mat3 = geometry::Mat3;

struct IG_GAttributes
{
	std::shared_ptr<IncidenceGraph::Attribute<Vec3>> vertex_position;
	std::shared_ptr<IncidenceGraph::Attribute<Scalar>> vertex_radius;

	std::shared_ptr<IncidenceGraph::Attribute<Vec3>> face_normal;

	std::shared_ptr<IncidenceGraph::Attribute<Dart>> vertex_contact_surface;
	std::shared_ptr<IncidenceGraph::Attribute<std::pair<Dart, Dart>>> halfedge_contact_surface_face;
	std::shared_ptr<IncidenceGraph::Attribute<std::pair<Mat3, Mat3>>> halfedge_frame;
	// std::shared_ptr<IncidenceGraph::Attribute<std::pair<Dart, Dart>>> halfedge_volume_connection;

	std::shared_ptr<IncidenceGraph::Attribute<std::vector<Dart>>> face_edge_dart;
	std::shared_ptr<IncidenceGraph::Attribute<std::vector<Dart>>> vertex_boundary_edge_dart;

	// std::shared_ptr<IncidenceGraph::Attribute<Dart>> vertex_up_dart;
	// std::shared_ptr<IncidenceGraph::Attribute<Vec3>> vertex_normal;
};

struct IG_M2Attributes
{
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_position;
	std::shared_ptr<CMap2::Attribute<std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge>>> dual_vertex_graph_branch;
	std::shared_ptr<CMap2::Attribute<IncidenceGraph::Vertex>> volume_igvertex;
	std::shared_ptr<CMap2::Attribute<Vec3>> volume_center;
	std::shared_ptr<CMap2::Attribute<Vec3>> edge_mid;
	std::shared_ptr<CMap2::Attribute<Dart>> halfedge_volume_connection;
	std::shared_ptr<CMap2::Attribute<CMap2*>> ortho_scaffold;
};

struct IG_M3Attributes
{
	std::shared_ptr<CMap3::Attribute<Vec3>> vertex_position;
	// std::shared_ptr<CMap3::Attribute<Graph::HalfEdge>> volume_graph_connection;

	DartMarker<CMap3>* extremity_faces;

	// std::shared_ptr<CMap3::Attribute<Mat3>> corner_frame;
	// std::shared_ptr<CMap3::Attribute<Mat3>> hex_frame;

	// std::shared_ptr<CMap3::Attribute<Scalar>> scaled_jacobian;
	// std::shared_ptr<CMap3::Attribute<Scalar>> jacobian;
	// std::shared_ptr<CMap3::Attribute<Scalar>> max_frobenius;
	// std::shared_ptr<CMap3::Attribute<Scalar>> mean_frobenius;

	// std::shared_ptr<CMap3::Attribute<Vec3>> color_scaled_jacobian;
	// std::shared_ptr<CMap3::Attribute<Vec3>> color_jacobian;
	// std::shared_ptr<CMap3::Attribute<Vec3>> color_max_frobenius;
	// std::shared_ptr<CMap3::Attribute<Vec3>> color_mean_frobenius;
};

std::tuple<IG_GAttributes, IG_M2Attributes, IG_M3Attributes> incidenceGraph_to_hex(IncidenceGraph& ig, CMap2& m2,
																				   CMap3& m3);

/*****************************************************************************/
/* data preparation                                                          */
/*****************************************************************************/

bool add_incidenceGraph_attributes(IncidenceGraph& ig, IG_GAttributes& igAttributes);
bool add_cmap2_attributes(CMap2& m2, IG_M2Attributes& m2Attribs);
bool add_cmap3_attributes_igh(CMap3& m3, IG_M3Attributes& m3Attribs);

bool compute_faces_geometry(const IncidenceGraph& ig, const IncidenceGraphData& incidenceGraph_data,
							IG_GAttributes& igAttributes);

/*****************************************************************************/
/* utils                                                                     */
/*****************************************************************************/

std::vector<IncidenceGraph::Vertex> get_branch_vertices(const IncidenceGraph& g, const Branch& b);
std::vector<IncidenceGraph::Edge> get_incident_leaflet_edges(const IncidenceGraph& ig, IncidenceGraph::Vertex v0,
															 IncidenceGraph::Edge e0);
IncidenceGraph::Edge get_shared_edge(const IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Vertex v1);

Dart add_chunk(CMap3& m3);
Dart add_plate(CMap3& m3);

void index_volume_cells_igh(CMap2& m, CMap2::Volume vol);
void sew_volumes_igh(CMap3& m, Dart d0, Dart d1);

Mat3 rmf_step(const Vec3& x0, const Vec3& x1, const Mat3& U0, const Vec3& t1);
Vec3 slerp_igh(Vec3 A, Vec3 B, Scalar coef, bool in);
Scalar angle_on_sphere_igh(Vec3 A, Vec3 B, Vec3 C);
Scalar edge_max_angle_igh(CMap2& m2, CMap2::Edge e, IG_M2Attributes& m2Attribs);
Scalar min_cut_angle_igh(CMap2& m2, CMap2::Vertex v0, CMap2::Vertex v1, IG_M2Attributes& m2Attribs);
Vec3 spherical_barycenter_igh(std::vector<Vec3>& points, uint32 iterations);

bool dijkstra_topo_igh(CMap2& m2, CMap2::Vertex v0, std::shared_ptr<CMap2::Attribute<Dart>> previous,
					   std::shared_ptr<CMap2::Attribute<uint32>> dist);
Dart convex_hull_around_vertex(const IncidenceGraph& g, IncidenceGraph::Vertex v, CMap2& m2, IG_M2Attributes& m2Attribs,
							   std::vector<Vec3>& Ppos, std::vector<IncidenceGraph::Edge>& Pev);
void dualize_volume(CMap2& m, CMap2::Volume vol, IG_M2Attributes& m2Attribs, const IncidenceGraph& ig,
					IG_GAttributes& igAttribs);
Dart remesh_igh(CMap2& m2, CMap2::Volume vol, IG_M2Attributes& m2Attribs);

/*****************************************************************************/
/* contact surfaces generation                                               */
/*****************************************************************************/

bool build_contact_surfaces(const IncidenceGraph& ig, IG_GAttributes& igAttribs,
							IncidenceGraphData& incidenceGraph_data, CMap2& m2, IG_M2Attributes& m2Attribs);
void build_contact_surface_1(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2, IG_M2Attributes& m2Attribs,
							 IncidenceGraph::Vertex v);
void build_contact_surface_2(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2, IG_M2Attributes& m2Attribs,
							 IncidenceGraph::Vertex v, uint32 nb_leaflets);
void build_contact_surface_orange(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
								  IG_M2Attributes& m2Attribs, IncidenceGraph::Vertex v, std::vector<Vec3>& Ppos,
								  std::vector<IncidenceGraph::Edge>& Pev);
void build_contact_surface_n(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2, IG_M2Attributes& m2Attribs,
							 IncidenceGraph::Vertex v, uint32 nb_leaflets, uint32 nb_edges);
bool build_contact_surface_ortho(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
								 IG_M2Attributes& m2Attribs, IncidenceGraph::Vertex v, uint32 nb_edges);

/*****************************************************************************/
/* frames initialization & propagation                                       */
/*****************************************************************************/

bool create_intersection_frames(const IncidenceGraph& ig, IG_GAttributes& igAttribs, const IncidenceGraphData& igData,
								CMap2& m2, IG_M2Attributes m2Attribs);
bool create_intersection_frame_n(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
								 IG_M2Attributes& m2Attribs, IncidenceGraph::Vertex v);
bool create_ef_frame(const IncidenceGraph& ig, IG_GAttributes& igAttribs, IncidenceGraph::Vertex v);
bool create_ff_frame(const IncidenceGraph& ig, IG_GAttributes& igAttribs, IncidenceGraph::Vertex v);
bool create_extremity_frame(const IncidenceGraph& ig, IG_GAttributes& igAttribs, IncidenceGraph::Vertex v);

bool propagate_frames(const IncidenceGraph& ig, IG_GAttributes& gAttribs, const IncidenceGraphData& igData, CMap2& m2);
void propagate_frame_n_1(const IncidenceGraph& ig, IG_GAttributes& igAttribs,
						 std::vector<IncidenceGraph::Vertex>& branch_vertices);
bool propagate_frame_n_n(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
						 std::vector<IncidenceGraph::Vertex>& branch_vertices);

/*****************************************************************************/
/* contact surfaces geometry                                                 */
/*****************************************************************************/

bool set_contact_surfaces_geometry(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
								   IG_M2Attributes& m2Attribs);

/*****************************************************************************/
/* volume mesh generation                                                    */
/*****************************************************************************/

bool build_volumes(IncidenceGraph& ig, IG_GAttributes& igAttribs, IncidenceGraphData& incidenceGraph_data, CMap2& m2,
				   IG_M2Attributes& m2Attribs, CMap3& m3, IG_M3Attributes& m3Attribs);
void insert_ortho_chunks(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2, IG_M2Attributes& m2Attribs,
						 CMap3& m3);
bool build_branch_section(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2, IG_M2Attributes& m2Attribs,
						  CMap3& m3, IG_M3Attributes& m3Attribs, IncidenceGraph::Edge e0);
bool build_leaflets(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2, IG_M2Attributes& m2Attribs, CMap3& m3,
					const IncidenceGraphData& incidenceGraph_data);
bool build_leaflet_plates(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap3& m3,
						  const std::vector<IncidenceGraph::Face>& leaflet);
bool sew_leaflet_fan_edge_plates(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap3& m3, IncidenceGraph::Edge e);
bool build_leaflet_boundary_edge_plate(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
									   IG_M2Attributes& m2Attribs, CMap3& m3, IncidenceGraph::Edge be);
bool sew_leaflet_boundary_vertex_corner_plates(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap3& m3,
											   IncidenceGraph::Vertex bv);
bool sew_leaflet_boundary_vertex_fan_plates(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap3& m3,
											IncidenceGraph::Vertex bv);
bool register_leaflet_boundary_edge_plates_into_efjuncture(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
														   IG_M2Attributes& m2Attribs, CMap3& m3,
														   IncidenceGraph::Edge be, IncidenceGraph::Vertex bv);
bool register_leaflet_boundary_edge_plates_into_ffjuncture(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
														   IG_M2Attributes& m2Attribs, CMap3& m3,
														   IncidenceGraph::Edge be, IncidenceGraph::Vertex bv);
bool register_leaflet_boundary_edge_plates_into_intersection(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
															 IG_M2Attributes& m2Attribs, CMap3& m3,
															 IncidenceGraph::Edge be, IncidenceGraph::Vertex bv);
bool sew_sections_igh(CMap2& m2, IG_M2Attributes& m2Attribs, CMap3& m3);
bool set_volumes_geometry_igh(IncidenceGraph& ig, IG_GAttributes& igAttribs, IncidenceGraphData& incidenceGraph_data,
							  CMap2& m2, IG_M2Attributes& m2Attribs, CMap3& m3, IG_M3Attributes& m3Attribs);

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_ALGOS_INCIDENCEGRAPH_TO_HEX_H_

// add_graph_attributes(g, gAttribs);
// add_cmap2_attributes(m2, m2Attribs);
// build_contact_surfaces(g, gAttribs, m2, m2Attribs);
// create_intersection_frames(g, gAttribs, m2, m2Attribs);
// propagate_frames(g, gAttribs, gData, m2);
// set_contact_surfaces_geometry(g, gAttribs, m2, m2Attribs);
// build_branch_sections(g, gAttribs, m2, m2Attribs, m3);
// sew_branch_sections(m2, m2Attribs, m3);

// add_cmap3_attributes(m3, m3Attribs);
// set_volumes_geometry(m2, m2Attribs, m3, m3Attribs);

// - analyse du graphe:
// 	+ détection des feuillets (+ configurations invalides)
// 		(config invalides = éventail && intersection branche surface variétée)
// 	- détection des branches avec leurs extrémités

// + calcul des normales des faces + moyennes de arêtes en coin

// - construction des échaffaudages
// 	- arêtes - valence 1/2: quad topo
// 	+ arete + coin de face || coin + coin: quad topo + géometrie
// 	- aretes valence >= 3 || arêtes + faces || >=3 coins de faces: partition de sphere -> topo + geometrie

// - préparation de la géometrie des feuillets
// - construction des repères à propager
// 	- depuis les embranchements complexes
// 	+ depuis les faces

// - propagation de la geometrie dans les échafaudages

// - construction des hexas
// 	- construction des troncons 4 hexs
// 	+ construction des palettes/plates (blocs d'hex sur les faces)

// - couture des hexs

// - plongement géometrique

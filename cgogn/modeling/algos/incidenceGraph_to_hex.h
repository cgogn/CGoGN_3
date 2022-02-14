// get_graph_data(g, gData);

#ifndef CGOGN_MODELING_ALGOS_INCIDENCEGRAPH_TO_HEX_H_
#define CGOGN_MODELING_ALGOS_INCIDENCEGRAPH_TO_HEX_H_

#include <cgogn/core/types/incidence_graph/incidence_graph_ops.h>
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

template <typename MESH, typename CELL>
class CellMarker;

namespace modeling
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;
using Mat3 = geometry::Mat3;


struct IGAttributes
{
	std::shared_ptr<IncidenceGraph::Attribute<Vec3>> vertex_position;
	std::shared_ptr<IncidenceGraph::Attribute<std::vector<uint8>>> incident_edges_dir;

	std::shared_ptr<IncidenceGraph::Attribute<Vec3>> face_normal;
	std::shared_ptr<IncidenceGraph::Attribute<Vec3>> face_center;
	std::shared_ptr<IncidenceGraph::Attribute<std::vector<Vec3>>> face_vertex_tangent;
	std::shared_ptr<IncidenceGraph::Attribute<Scalar>> vertex_radius;
	// std::shared_ptr<IncidenceGraph::Attribute<Dart>> vertex_contact_surface;
	std::shared_ptr<IncidenceGraph::Attribute<std::vector<Dart>>> vertex_boundary_edge_dart;
	std::shared_ptr<IncidenceGraph::Attribute<Dart>> vertex_contact_surface;

	std::shared_ptr<IncidenceGraph::Attribute<std::pair<Dart, Dart>>> halfedge_volume_connection;
	std::shared_ptr<IncidenceGraph::Attribute<std::vector<Dart>>> face_edge_dart;
	std::shared_ptr<IncidenceGraph::Attribute<std::pair<Dart, Dart>>> halfedge_contact_surface_face;
	std::shared_ptr<IncidenceGraph::Attribute<std::pair<Mat3, Mat3>>> halfedge_frame;

	std::shared_ptr<IncidenceGraph::Attribute<Dart>> vertex_up_dart;
	std::shared_ptr<IncidenceGraph::Attribute<Vec3>> vertex_normal;

};

struct M2Attributes
{
	std::shared_ptr<CMap2::Attribute<Vec3>> vertex_position;
	std::shared_ptr<CMap2::Attribute<std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge>>> dual_vertex_graph_branch;
	std::shared_ptr<CMap2::Attribute<IncidenceGraph::Vertex>> volume_igvertex;
	std::shared_ptr<CMap2::Attribute<Vec3>> volume_center;
	std::shared_ptr<CMap2::Attribute<Vec3>> edge_mid;
	std::shared_ptr<CMap2::Attribute<Dart>> halfedge_volume_connection;
	// std::shared_ptr<CMap2::Attribute<CMap2*>> ortho_scaffold;
};

struct M3Attributes
{
	std::shared_ptr<CMap3::Attribute<Vec3>> vertex_position;
	// std::shared_ptr<CMap3::Attribute<Graph::HalfEdge>> volume_graph_connection;

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


struct IncidenceGraphData
{
	std::vector<std::pair<IncidenceGraph::Edge, IncidenceGraph::Edge>> branches;
	std::vector<IncidenceGraph::Vertex> efjunctures;
	std::vector<IncidenceGraph::Vertex> ffjunctures;
	std::vector<IncidenceGraph::Vertex> intersections;
	std::vector<IncidenceGraph::Face> leaflets;
};

// struct VertexData
// {
// 	uint32 degree;
// 	uint32 isolatedEdge;
// };

// inline uint32 pseudoDegree(IncidenceGraph& ig, IncidenceGraph::Vertex v)
// {
// 	uint32 degree = 0;

// 	foreach_incident_edge(ig, v, [&](IncidenceGraph::Edge e) -> bool {
// 		uint32 count = 0;
// 		foreach_incident_face(ig, e, [&](IncidenceGraph::Face f) -> bool {
// 			++count;
// 			return true;
// 		});
		
// 		switch(count)
// 		{
// 			case 0: 
// 				++degree;
// 				break;
// 			case 1:
// 				degree += 0.5;
// 				break;
// 			case 2:
// 				break;
// 			default:
// 				degree = -1;
// 				break;
// 		}

// 		return (degree != -1);
// 	});

// 	return degree;
// }

std::tuple<IGAttributes, M2Attributes, M3Attributes> incidenceGraph_to_hex(IncidenceGraph& ig, CMap2& m2, CMap3& m3);

/*****************************************************************************/
/* data preparation                                                          */
/*****************************************************************************/

bool get_incidenceGraph_data(const IncidenceGraph& ig, IncidenceGraphData& incidenceGraph_data);
bool add_incidenceGraph_attributes(IncidenceGraph& ig, IGAttributes& igAttributes);
bool add_cmap2_attributes(CMap2& m2, M2Attributes& m2Attribs);
bool compute_faces_geometry(const IncidenceGraph& ig, const IncidenceGraphData& incidenceGraph_data, IGAttributes& igAttributes);
bool add_cmap3_attributes_igh(CMap3& m3, M3Attributes& m3Attribs);

/*****************************************************************************/
/* utils                                                                     */
/*****************************************************************************/

std::pair<IncidenceGraph::Edge, IncidenceGraph::Edge> find_branch_extremities(const IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Edge e0, CellMarker<IncidenceGraph, IncidenceGraph::Edge>& cm);
std::vector<IncidenceGraph::Face> get_incident_leaflets(const IncidenceGraph& ig, IncidenceGraph::Vertex v0);
std::vector<IncidenceGraph::Face> get_leaflet_faces(const IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Face f0);
std::vector<IncidenceGraph::Vertex> get_branch_vertices(const IncidenceGraph& g, IncidenceGraph::Edge e0);
bool contains_vertex(const IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Face f0);
std::vector<IncidenceGraph::Face> get_incident_leaflet(const IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Face f0);
std::vector<IncidenceGraph::Edge> get_incident_leaflet_edges(const IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Edge e0);
std::vector<IncidenceGraph::Edge> get_leaflet_boundary_edges(const IncidenceGraph& ig, std::vector<IncidenceGraph::Face>& leaflet);
std::vector<IncidenceGraph::Edge> get_leaflet_boundary_edges(const IncidenceGraph& ig, IncidenceGraph::Face leaflet);
std::vector<IncidenceGraph::Vertex> get_leaflet_boundary_vertices(const IncidenceGraph& ig, std::vector<IncidenceGraph::Face>& leaflet);
std::vector<IncidenceGraph::Face> get_leaflets(const IncidenceGraph& ig);
Dart add_chunk(CMap3& m3);
Dart add_plate(CMap3& m3);

void rmf(const IncidenceGraph& ig, IGAttributes& igAttribs, std::vector<IncidenceGraph::Vertex> branch_vertices);
Mat3 rmf_step(Vec3 x0, Vec3 x1, Mat3 U0, Vec3 t1);

void index_volume_cells_igh(CMap2& m, CMap2::Volume vol);
void sew_volumes_igh(CMap3& m, Dart d0, Dart d1);
bool dijkstra_topo_igh(CMap2& m2, CMap2::Vertex v0, std::shared_ptr<CMap2::Attribute<Dart>> previous,
				   std::shared_ptr<CMap2::Attribute<uint32>> dist);
Dart convex_hull_around_vertex(const IncidenceGraph& g, IncidenceGraph::Vertex v, CMap2& m2, M2Attributes& m2Attribs,
							   std::vector<Vec3>& Ppos, std::vector<IncidenceGraph::Edge>& Pev);
void dualize_volume(CMap2& m, CMap2::Volume vol, M2Attributes& m2Attribs, const IncidenceGraph& ig, IGAttributes& igAttribs);
Dart remesh_igh(CMap2& m2, CMap2::Volume vol, M2Attributes& m2Attribs);
Vec3 slerp_igh(Vec3 A, Vec3 B, Scalar coef, bool in);
Scalar angle_on_sphere_igh(Vec3 A, Vec3 B, Vec3 C);
Scalar edge_max_angle_igh(CMap2& m2, CMap2::Edge e, M2Attributes& m2Attribs);
Scalar min_cut_angle_igh(CMap2& m2, CMap2::Vertex v0, CMap2::Vertex v1, M2Attributes& m2Attribs);
Vec3 spherical_barycenter_igh(std::vector<Vec3>& points, uint32 iterations);

uint32 get_incident_edge_id(const IncidenceGraph& ig, IncidenceGraph::Vertex v, IncidenceGraph::Edge e);
uint32 get_incident_edge_id(const IncidenceGraph& ig, IncidenceGraph::Face f, IncidenceGraph::Edge e);
uint32 get_incident_vertex_id(const IncidenceGraph& ig, IncidenceGraph::Face f, IncidenceGraph::Vertex v);
IncidenceGraph::Edge get_shared_edge(const IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Vertex v1);


/*****************************************************************************/
/* contact surfaces generation                                               */
/*****************************************************************************/


bool build_contact_surfaces(const IncidenceGraph& ig, IGAttributes& igAttribs, IncidenceGraphData& incidenceGraph_data, CMap2& m2, M2Attributes& m2Attribs);
void build_contact_surface_1(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs,
							 IncidenceGraph::Vertex v);
void build_contact_surface_2(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs,
							 IncidenceGraph::Vertex v);
void build_contact_surface_n(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs, IncidenceGraph::Vertex v);
void build_contact_surface_orange(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs,
								  IncidenceGraph::Vertex v, std::vector<Vec3>& Ppos, std::vector<IncidenceGraph::Edge>& Pev);
// bool build_contact_surface_ortho(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs,
// 								 Graph::Vertex v);

/*****************************************************************************/
/* frames initialization & propagation                                       */
/*****************************************************************************/

bool create_intersection_frames(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes m2Attribs);
bool create_intersection_frame_n(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs,
								 IncidenceGraph::Vertex v);
bool create_ef_frame(const IncidenceGraph& ig, IGAttributes& igAttribs, IncidenceGraph::Vertex v);
bool create_extremity_frame(const IncidenceGraph& ig, IGAttributes& igAttribs, IncidenceGraph::Vertex v);

bool propagate_frames(const IncidenceGraph& ig, IGAttributes& gAttribs, const IncidenceGraphData& igData, CMap2& m2);
void propagate_frame_n_1(const IncidenceGraph& ig, IGAttributes& igAttribs,  std::vector<IncidenceGraph::Vertex>& branch_vertices);
bool propagate_frame_n_n(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2,  std::vector<IncidenceGraph::Vertex>& branch_vertices);


/*****************************************************************************/
/* contact surfaces geometry                                                 */
/*****************************************************************************/
bool prepare_leaflets_geometry(const IncidenceGraph& ig, IGAttributes& igAttribs, IncidenceGraphData& incidenceGraph_data, CMap3& m3);
bool set_contact_surfaces_geometry(const IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs);

/*****************************************************************************/
/* volume mesh generation                                                    */
/*****************************************************************************/
bool build_volumes(IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3);
bool build_branch_section(IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3, IncidenceGraph::Edge e0);
bool build_leaflet(IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3, IncidenceGraph::Face f0);
bool build_leaflet_edge(IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3, IncidenceGraph::Face f0);
bool sew_sections_igh(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3);
bool set_volumes_geometry_igh(IncidenceGraph& ig, IGAttributes& igAttribs, CMap2& m2, M2Attributes& m2Attribs, CMap3& m3, M3Attributes& m3Attribs);
// bool set_volumes_geometry_igh(CMap2& m2, M2Attributes& m2Attribs, CMap3& m3, M3Attributes& m3Attribs);


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
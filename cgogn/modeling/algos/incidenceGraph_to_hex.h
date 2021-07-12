// get_graph_data(g, gData);

#ifndef CGOGN_MODELING_ALGOS_INCIDENCEGRAPH_TO_HEX_H_
#define CGOGN_MODELING_ALGOS_INCIDENCEGRAPH_TO_HEX_H_

#include <cgogn/core/types/incidence_graph/incidence_graph_ops.h>

namespace cgogn
{

template <typename MESH, typename CELL>
class CellMarker;

namespace modeling
{

struct IncidenceGraphData
{
	std::vector<std::pair<IncidenceGraph::Edge, IncidenceGraph::Edge>> branches;
	std::vector<Graph::Vertex> efjunctures;
	std::vector<Graph::Vertex> ffjunctures;
	std::vector<Graph::Vertex> intersections;
	std::vector<IncidenceGraph::Face> leaflets;
};

bool get_incidenceGraph_data(const MESH& ig, IncidenceGraphData& incidenceGraph_data);


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

// - construction des repères à propager
// 	- depuis les embranchements complexes
// 	+ depuis les faces

// - propagation de la geometrie dans les échafaudages

// - construction des hexas
// 	- construction des troncons 4 hexs
// 	+ construction des palettes/plates (blocs d'hex sur les faces)

// - couture des hexs 

// - plongement géometrique
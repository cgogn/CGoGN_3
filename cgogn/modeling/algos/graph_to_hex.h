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
	std::shared_ptr<CMap2::Attribute<Vec3>> volume_center;
	std::shared_ptr<CMap2::Attribute<Vec3>> edge_mid;
	std::shared_ptr<CMap2::Attribute<Dart>> halfedge_volume_connection;
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
bool dijkstra_topo(CMap2& m2, std::shared_ptr<CMap2::Attribute<Dart>> previous, std::shared_ptr<CMap2::Attribute<uint>> dist);
Dart convex_hull(CMap2& m2, const cgogn::io::SurfaceImportData& surface_data);
Dart remesh(CMap2& m, CMap2::Volume vol, M2Attributes& m2Attribs);
Vec3 slerp(Vec3 A, Vec3 B, Scalar coef, bool in);
Scalar angle_on_sphere(Vec3 A, Vec3 B, Vec3 C);
Scalar edge_max_angle(CMap2& m2, CMap2::Edge e, M2Attributes& m2Attribs);
Scalar min_cut_angle(CMap2& m2, CMap2::Vertex v0, CMap2::Vertex v1, M2Attributes& m2Attribs);

/*****************************************************************************/
/* data preparation                                                          */
/*****************************************************************************/

bool get_graph_data(const Graph& g, GData& gData);
bool add_graph_attributes(Graph& g, GAttributes& gAttribs);
bool add_cmap2_attributes(CMap2& m2, M2Attributes& m2Attribs);

/*****************************************************************************/
/* contact surfaces generation                                               */
/*****************************************************************************/

bool build_contact_surfaces(const Graph& g, GAttributes& gAttribs, CMap2& m2, M2Attributes& m2Attribs);
void build_contact_surface_1(const Graph& g, GAttributes& gAttribs, CMap2& m2, Graph::Vertex v);
void build_contact_surface_2(const Graph& g, GAttributes& gAttribs, CMap2& m2, Graph::Vertex v);
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

// Dart convex_quad(CMap2& m2, M2Attributes& m2Attribs, Dart f){
//         Dart res;
//         Vec3 A = value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(f));
//         Vec3 B = value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(m2.phi1(f)));
//         Vec3 C = value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(m2.phi<11>(f)));
//         Vec3 D = value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(m2.phi_1(f)));
//         Vec3 AB = B - A;
//         Vec3 AC = C - A;
//         Vec3 AD = D - A;
//         Vec3 N0 = AC.cross(AD);
//         Vec3 N1 = AB.cross(AC);
//         Vec3 N = N0.cross(N1).normalized();

//         if (N.dot(AC) < 0){
//             res = m2.phi1(f);
//         }
//         else{
//             res = f;
//         }
//         return res;
//     }

// Vec3 mean_dir(Vec3 center, Scalar radius, Vec3 point, std::vector<Vec3> points){
//     using Quat = Eigen::Quaterniond;
//     uint32 valence = points.size();

//     std::vector<Vec3> directions;
//     for (Vec3 p : points){
//         directions.push_back((p - center).normalized());
//     }
//     Vec3 avg_dir = (point - center).normalized();

//     std::vector<Quat> rotations;
//     rotations.reserve(valence);

// //        for (uint i = 0; i < valence; i++){
//         for (Vec3 dir : directions){
//             Quat q = Quat::FromTwoVectors(avg_dir, dir);
//             q.normalize();
//             rotations.push_back(q);
//         }

//         Eigen::MatrixXd m(4, valence);
//         for (uint32 j = 0; j < valence; ++j){
//             const Quat& q = rotations[j];
//             m.col(j) = Eigen::Vector4d(q.w(), q.x(), q.y(), q.z());
//         }

//         Eigen::MatrixXd mm = m * m.transpose();
//         Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mm);
//         Eigen::Vector4d r = es.eigenvectors().col(3);

//         Quat mean_rot(r[0], r[1], r[2], r[3]);
//         mean_rot.normalize();
//         avg_dir = mean_rot._transformVector(avg_dir);
//         rotations.clear();
// //        }

//     return avg_dir * radius + center;
// }

// bool complete_intersection_n(
// 	const Graph& g,
// 	GAttributes& gAttribs,
//     CMap2& m2,
// 	M2Attributes& m2Attribs,
// 	Graph::Vertex gv
// )
// {
//     using geometry::intersection_ray_triangle;

//     Vec3 center = value<Vec3>(g, gAttribs.vertex_position, gv);
//     Scalar radius = value<Scalar>(g, gAttribs.vertex_radius, gv);
//     Dart cc = value<Dart>(g, gAttribs.m2_CC, gv);

//     std::vector<Dart> Fd = {cc, m2.phi<2111>(cc), m2.phi<12>(cc)};
//     std::vector<Vec3> Ppos;
//     Ppos.reserve(3);
//     std::vector<Dart> Pdart;
//     Pdart.reserve(3);
//     g.foreach_dart_of_orbit(gv, [&](Dart d) -> bool
//     {
//         Vec3 p = value<Vec3>(g, gAttribs.vertex_position, Graph::Vertex(g.alpha0(d)));
//         p = project_on_sphere(p, radius, center);
//         Ppos.push_back(p);
//         Pdart.push_back(d);
//         return true;
//     });

//     Vec3 V = (Ppos[1] - Ppos[0]).cross(Ppos[2] - Ppos[0]).normalized();
//     std::vector<Vec3> Q = {center + V * radius, center - V * radius};
//     std::vector<Vec3> M = {center + (Ppos[1] - Ppos[0]).normalized().cross(V) * radius,
//         center + (Ppos[2] - Ppos[1]).normalized().cross(V) * radius,
//         center + (Ppos[0] - Ppos[2]).normalized().cross(V) * radius};

//     value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(Fd[0])) = M[0];
//     value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(Fd[1])) = M[1];
//     value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(Fd[2])) = M[2];
//     value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(m2.phi1(Fd[0]))) = Q[0];
//     value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(m2.phi_1(Fd[0]))) = Q[1];

//     value<Vec3>(m2, m2Attribs.face_point, M2Face(Fd[0])) = Ppos[0];
//     value<Vec3>(m2, m2Attribs.face_point, M2Face(Fd[1])) = Ppos[1];
//     value<Vec3>(m2, m2Attribs.face_point, M2Face(Fd[2])) = Ppos[2];

//     value<Dart>(m2, m2Attribs.face_branch, M2Face(Fd[0])) = Pdart[0];
//     value<Dart>(m2, m2Attribs.face_branch, M2Face(Fd[1])) = Pdart[1];
//     value<Dart>(m2, m2Attribs.face_branch, M2Face(Fd[2])) = Pdart[2];

//     for (uint i = 3; i < Ppos.size(); i++)
//     {
//         Vec3 P0 = Ppos[i];
//         Dart F0 = Pdart[i];
//         M2Face face2cut;

//         std::vector<Vec3> Quadp;
//         std::vector<Dart> Quadd;

//         bool face_found = false;

//         foreach_cell(m2, [&](M2Face f) -> bool
//         {
//             Dart fd = convex_quad(m2, m2Attribs, f.dart);
//             Quadd = {fd, m2.phi1(fd),
//                         m2.phi<11>(fd),
//                         m2.phi_1(fd)};
//             Quadp = {value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(Quadd[0])),
//                     value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(Quadd[1])),
//                     value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(Quadd[2])),
//                     value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(Quadd[3]))};
//             if ((face_found = (intersection_ray_triangle(center, P0 - center, Quadp[0], Quadp[1], Quadp[2])
//                                 || intersection_ray_triangle(center, P0 - center, Quadp[0], Quadp[2], Quadp[3])))){
//                 face2cut = f;
//                 return false;
//             }
//             return true;
//         });
//         if (!face_found)
//         {
//             std::cout << "error complete_intersection_n: no face found" << std::endl;
//             return false;
//         }

//         Vec3 P1 = value<Vec3>(m2, m2Attribs.face_point, M2Face(face2cut));
//         Dart F1 = value<Dart>(m2, m2Attribs.face_branch, M2Face(face2cut));

//         Dart cut0, cut1;
//         Vec3 AC = (Quadp[2] - Quadp[0]).normalized();
//         Vec3 BD = (Quadp[3] - Quadp[1]).normalized();
//         Vec3 P0P1 = (P1 - P0).normalized();
//         if (abs(AC.dot(P0P1)) < abs(BD.dot(P0P1))){
//             cut0 = Quadd[0]; cut1 = Quadd[2];
//         } else {
//             cut0 = Quadd[1]; cut1 = Quadd[3];
//         }
//         cut_face(m2, M2Vertex(cut0), M2Vertex(cut1));
//         M2Vertex v = cut_edge(m2, M2Edge(m2.phi_1(cut0)));
//         value<Vec3>(m2, m2Attribs.vertex_position, v) = project_on_sphere((P0 + P1) * Scalar(0.5), radius, center);

//         Dart newFace, oldFace;
//         Vec3 out0 =  value<Vec3>(m2, m2Attribs.vertex_position, v) - center;
//         // Vec3 out1 = ((P0 - value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(cut0))).normalized().cross((P1 -
//         m2_attribs.position[cut0]).normalized()));

//         Vec3 out1 = ((P0 - value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(cut0))).normalized().cross((P1 -
//         value<Vec3>(m2, m2Attribs.vertex_position, M2Vertex(cut0))).normalized())); if (out1.dot(out0) >= 0)
//         {
//             newFace = cut0; oldFace = cut1;
//         }
//         else
//         {
//             newFace = cut1; oldFace = cut0;
//         }

//         value<Vec3>(m2, m2Attribs.face_point, M2Face(newFace)) = P0;
//         value<Vec3>(m2, m2Attribs.face_point, M2Face(oldFace)) = P1;
//         value<Dart>(m2, m2Attribs.face_branch, M2Face(newFace)) = F0;
//         value<Dart>(m2, m2Attribs.face_branch, M2Face(oldFace)) = F1;

//         foreach_incident_vertex(m2, M2Face(newFace), [&](M2Vertex m2v) -> bool
//         {
//             std::vector<Vec3> points;
//             foreach_incident_face(m2, m2v, [&](M2Face m2f) -> bool
//             {
//                 points.push_back(value<Vec3>(m2, m2Attribs.face_point, m2f));
//                 return true;
//             });
//             value<Vec3>(m2, m2Attribs.vertex_position, m2v) = mean_dir(center, radius, value<Vec3>(m2,
//             m2Attribs.vertex_position, m2v), points); return true;
//         });
//     }

//     foreach_incident_face(m2, M2Volume(cc), [&](M2Face m2f) -> bool
//     {
//         Dart f = value<Dart>(m2, m2Attribs.face_branch, m2f);
//         value<Dart>(g, gAttribs.m2_interface, GDart(f)) = m2f.dart;
//         return true;
//     });

//     return true;
// }

// bool complete_intersections(
// 	const Graph& g,
// 	GAttributes& gAttribs,
// 	CMap2& m2,
// 	M2Attributes& m2Attribs)
// {
// 	foreach_cell(g, [&] (Graph::Vertex v) -> bool
// 	{
// 		switch (degree(g, v))
// 		{
// 		case 1:
// 			break;
// 		case 2:
// 			break;
// 		case 3:
// 			complete_intersection_3(g, gAttribs, m2, m2Attribs, v);
// 			break;
// 		default:
// 			// complete_intersection_n(g, gAttribs, m2, m2Attribs, v);
// 			break;
// 		}
// 		return true;
// 	});

// 	return true;
// }

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_GRAPH_TO_HEX_H_

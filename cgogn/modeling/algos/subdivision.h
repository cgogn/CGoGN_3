/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_GEOMETRY_ALGOS_SUBDIVISION_H_
#define CGOGN_GEOMETRY_ALGOS_SUBDIVISION_H_

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/algos/angle.h>

#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_ops/volume.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/mr_cmap_infos.h>

namespace cgogn
{

namespace modeling
{

using Vec3 = geometry::Vec3;

///////////
// CMap2 //
///////////

void hexagon_to_triangles(CMap2& m, CMap2::Face f)
{
	cgogn_message_assert(codegree(m, f) == 6, "hexagon_to_triangles: given face should have 6 edges");
	Dart d0 = phi1(m,f.dart);
	Dart d1 = phi<11>(m,d0);
	cut_face(m, CMap2::Vertex(d0), CMap2::Vertex(d1));
	Dart d2 = phi<11>(m,d1);
	cut_face(m, CMap2::Vertex(d1), CMap2::Vertex(d2));
	Dart d3 = phi<11>(m,d2);
	cut_face(m, CMap2::Vertex(d2), CMap2::Vertex(d3));
}

//////////////
// MESHVIEW //
//////////////

template <typename MESH,
		  typename std::enable_if<is_mesh_view<MESH>::value>::type* = nullptr>
void
hexagon_to_triangles(MESH& m, typename mesh_traits<MESH>::Face f)
{
	hexagon_to_triangles(m.mesh(), f);
}

/////////////
// GENERIC //
/////////////

template <typename MESH>
void subdivide(MESH& m, typename mesh_traits<MESH>::template Attribute<Vec3>* vertex_position)
{
	using Vertex = typename cgogn::mesh_traits<MESH>::Vertex;
	using Edge = typename cgogn::mesh_traits<MESH>::Edge;
	using Face = typename cgogn::mesh_traits<MESH>::Face;

	CellCache<MESH> cache(m);
	cache.template build<Edge>();
	cache.template build<Face>();

	foreach_cell(cache, [&] (Edge e) -> bool
	{
		std::vector<Vertex> vertices = incident_vertices(m, e);
		Vertex v = cut_edge(m, e);
		value<Vec3>(m, vertex_position, v) =
			0.5 * (value<Vec3>(m, vertex_position, vertices[0]) + value<Vec3>(m, vertex_position, vertices[1]));
		return true;
	});

	foreach_cell(cache, [&] (Face f) -> bool
	{
		if (codegree(m, f) == 6)
			hexagon_to_triangles(m, f);
		return true;
	});
}

/* -------------------------- BUTTERFLY VOLUME MASKS -------------------------- */

// Remplit les vecteurs de p/q-points pour le volume incident à d
// Si le masque est mal définit (proche du bord) on ne garde que les p-points
// et le vecteur de q-points est vide après l'appel
template <typename MESH>
void volumePointMask(const MESH& m,Dart d, std::vector<Dart> &p_point, std::vector<Dart> &q_point)
{
    // p-points : ajouter les sommets de Face(d) + Face(phi<2112>(m,d))
    Dart b = d;
    Dart t = b;
    do {
        p_point.push_back(t);
        t = phi1(m,t);
    } while (t != b);
    b = phi<2112>(m,d);
    t = b;
    do {
        p_point.push_back(t);
        t = phi1(m,t);
    } while (t != b);

    // q-points : pour chaque face f du volume,
    // ajouter Face(Inherit::phi<32112>(f)) (face opposée du volume adjacent)
    // si le volume adjacent n'existe pas (bord) alors le masque est mal défini,
    // donc on ignore les q-points
    Dart f = d;
    unsigned i = 0;
    do {
        // S'il n'y a pas de q-points pour cette face,
        // on vide le vecteur de q-points et on quitte la recherche de q-points
        if (is_boundary(m,phi3(m,f))) {
            q_point.clear();
            break;
        }
        // Sinon, récupération des brins de la face opposée du volume adjacent
        b = phi<32112>(m,f);
        t = b;
        do {
            q_point.push_back(t);
            t = phi1(m,t);
        } while (t != b);
        // face suivante sur l'hexaèdre
        f = phi2(m,(i%2==0 ? phi1(m,f) : phi_1(m,f)));
        ++i;
    } while (f != d);
}

// Remplit les vecteurs de p/q/r/s/t-points pour la face incidente à d
// Si le masque est mal définit (proche du bord) on ne garde que les p-points
// et les autres vecteurs sont vides après l'appel
template <typename MESH>
void facePointMask(const MESH& m,Dart d, std::vector<Dart> &p_point, std::vector<Dart> &q_point,
                   std::vector<Dart> &r_point, std::vector<Dart> &s_point, std::vector<Dart> &t_point)
{
    // p-points : ajouter les sommets de Face(d)
    Dart t = d;
    foreach_incident_vertex(m,typename MESH::Face(d),[&](typename MESH::Vertex dd)->bool{
        p_point.push_back(dd.dart);
		return true;
    });

    // q-points : ajouter les sommets de Face(Inherit::phi<2112>(d)) + Face(Inherit::phi<32112>(d))
    Dart b = phi<2112>(m,d);
    t = b;
    do {
        q_point.push_back(t);
        t = phi1(m,t);
    } while (t != b);
    b = phi<32112>(m,d);
    t = b;
    do {
        q_point.push_back(t);
        t = phi1(m,t);
    } while (t != b);

    // r/t-points
    // Question : si un volume adjacent au volume sur ou sous la face n'existe pas, doit-on chercher
    // les r/s/t-points par d'autres chemins ? (nombre de config énorme...)

    // Pour chaque brin de p-point
    for (Dart p : p_point) {
     // Est-ce que les volumes adjacents aux volumes au dessus et au dessous existent ou est-ce en surface
        bool    face_phi2_is_surface = is_incident_to_boundary(m,typename MESH::Face(phi2(m,p))),
                face_phi32_is_surface = is_incident_to_boundary(m,typename MESH::Face(phi<32>(m,p)));

        Dart    phi323211 = phi<323211>(m,p),
                phi3232111 = phi<3232111>(m,p),
                phi23211 = phi<23211>(m,p),
                phi232111 = phi<232111>(m,p);

        // Si Face(phi<2>(p)) est sur la surface et Face(phi<32>(p)) aussi
        // alors il n'y a pas de r/t-points
        if (face_phi2_is_surface && face_phi32_is_surface) {
            q_point.clear();
            r_point.clear();
            s_point.clear();
            t_point.clear();
            break;
        }

        // Si Face(phi<2>(p)) est sur la surface et pas Face(phi<32>(p))
        // alors il y a deux r-points : Edge(phi<323211>(p))
        if (face_phi2_is_surface && !face_phi32_is_surface) {
            r_point.push_back(phi323211);
            r_point.push_back(phi3232111);
            continue;
        }

        // Si Face(phi<2>(p)) n'est pas sur la surface mais que Face(phi<32>(p)) l'est
        // alors il y a deux r-points : Edge(phi<23211>(p))
        if (!face_phi2_is_surface && face_phi32_is_surface) {
            r_point.push_back(phi23211);
            r_point.push_back(phi232111);
            continue;
        }

        // Ni Face(phi<2>(p)) ni Face(phi<32>(p)) ne sont sur la surface
        // alors soit Edge(phi<23211>(p)) == Edge(phi<323211>(p)) (topologie régulière)
        // soit Edge(phi<23211>(p)) != Edge(phi<323211>(p)) (arete extraordinaire)

        // Si Edge(phi<23211>(p)) != Edge(phi<323211>(p)) (arete extraordinaire)
        // alors il y a 4 t-points : Edge(phi<23211>(p)) + Edge(phi<323211>(p))
        if (index_of(m,typename MESH::Edge(phi23211)) != index_of(m,typename MESH::Edge(phi323211))) {
            t_point.push_back(phi23211);
            t_point.push_back(phi232111);
            t_point.push_back(phi323211);
            t_point.push_back(phi3232111);
        }
        // Sinon Edge(phi<23211>(p)) == Edge(phi<323211>(p))
        // et donc il y a 2 r-points
        else {
            r_point.push_back(phi23211);
            r_point.push_back(phi232111);
        }
    }

    // s-points : pour chaque q, ajouter Edge(phi<23211>(q))
    for (Dart q : q_point) {
        // S'il n'y a pas de s-points
        if (is_boundary(m,phi<23>(m,q))) {
            q_point.clear();
            r_point.clear();
            s_point.clear();
            t_point.clear();
            break;
        }
        s_point.push_back(phi<23211>(m,q));
        s_point.push_back(phi<232111>(m,q));
    }
}

// Remplit les vecteurs de p/q/r/s-points pour l'arete incidente à d
// Si le masque est mal définit (proche du bord) on ne garde que les p-points
// et les autres vecteurs sont vides après l'appel
// Attention, un appel récursif est fait (un seul)
template<typename MESH>
void edgePointMask(const MESH& m,Dart d, std::vector<Dart> &p_point, std::vector<Dart> &q_point,
                   std::vector<Dart> &r_point, std::vector<Dart> &s_point)
{
    if (is_incident_to_boundary(m,typename MESH::Edge(d))) {
        std::cerr << "Fatal error : edgePointMask called on dart which is surface" << std::endl;
        exit(EXIT_FAILURE);
    }

    // pour savoir si c'est le premier appel ou le rappel pour la deuxieme partie des r-points
    bool first_call = p_point.size() == 0;

    // p-points : ajouter les extremités de l'arete (sauf si c'est le deuxième appel
    if (first_call) {
        p_point.push_back(d);
        p_point.push_back(phi1(m,d));
    }

    // Pour toutes les faces incidentes à l'arete (au moins 4)
    // (arret garanti car l'arete n'est pas en surface)
    Dart t = d;
    do {
        // q-points : ajouter l'extremité de l'arete opposée (autre extremité : appel récursif)
        q_point.push_back(phi<11>(m,t));

        // r/s-points du coté sens direct
        // Question : s'il n'y a pas de volume adjacent à celui de la face courante,
        // doit-on chercher le point r/s par d'autres chemins (nombre de config!)

        bool face_phi12_is_surface = is_incident_to_boundary(m,typename MESH::Face(phi<12>(m,t))),
             face_phi132_is_surface = is_incident_to_boundary(m,typename MESH::Face(phi<132>(m,t)));

        Dart phi1323211 = phi<1323211>(m,t),
             phi1232_1 = phi_1(m,phi<1232>(m,t));

        // Si Face(phi<12>(t)) et Face(phi<132>(t)) sont sur la surface
        // alors pas de r/s-points coté direct et on peux donc arreter la recherche le points
        if (face_phi12_is_surface && face_phi132_is_surface) {
            q_point.clear();
            r_point.clear();
            s_point.clear();
            break;
        }
        // Au moins l'une des face n'est pas sur la surface
        // Si ni Face(phi<12>(t)) ni Face(phi<132>(t)) ne sont sur la surface
        // alors soit Vertex(phi<1323211>(t)) == Vertex(phi<1232-1>(t)) (topologie régulière)
        // soit Vertex(phi<1323211>(t)) != Vertex(phi<1232-1>(t)) (arete extraordinaire)
        if (!face_phi12_is_surface && !face_phi132_is_surface) {
            // Si Vertex(phi<1323211>(t)) != Vertex(phi<1232-1>(t)) (arete extraordinaire)
            // alors 2 s-points Vertex(phi<1323211>(t)) et Vertex(phi<1232-1>(t))
            if (index_of(m, typename MESH::Vertex(phi1323211)) != index_of(m, typename MESH::Vertex(phi1232_1))) {
                s_point.push_back(phi1323211);
                s_point.push_back(phi1232_1);
            }
            // Sinon Vertex(phi<1323211>(t)) == Vertex(phi<1232-1>(t)) (topologie régulière)
            // alors 1 r-point
            else {
                r_point.push_back(phi1323211);
            }
        }
        // sinon si Face(phi<12>(t)) est sur la surface (et pas Face(phi<132>(t)))
        // alors 1 r-point : Vertex(phi<1323211>(t))
        else if (face_phi12_is_surface) {
            r_point.push_back(phi1323211);
        }
        // sinon (si Face(phi<132>(t)) est sur la surface et pas Face(phi<12>(t))))
        // alors 1 r-point : Vertex(phi<1232-1>(t))
        else {
            r_point.push_back(phi1232_1);
        }

        // Face incidente suivante
        t = phi<23>(m,t);
    } while (t != d);

    // Sinon, faire la meme chose dans l'autre sens (si ce n'est pas déjà fait)
    // pour avoir les r/s-points de l'autre coté de l'arete
    if (q_point.size() != 0 && first_call) {
        edgePointMask(m,phi2(m,d), p_point, q_point, r_point, s_point);
    }
}

/* ------------------------- BUTTERFLY SURFACE MASKS ------------------------- */

// retourne le brin sur la même arete que d (en surface) mais sur la face adjacente en surface toujours
template<typename MESH>
Dart surfaceAdjacentFace(const MESH& m,Dart d)
{
    if (!is_boundary(m,phi3(m,d))) {
        std::cerr << "Fatal error : surfaceAdjacentFace called on dart which is not surface" << std::endl;
        exit(EXIT_FAILURE);
    }
    Dart t = phi2(m,d);
    // Recherche de la face adjacente (arrêt garanti car d est en surface)
    while (!is_boundary(m,phi3(m,t))) {
        t = phi<32>(m,t);
    }
    return t;
}

// remplit les vecteurs de p/q-points pour la face incidente à d à la surface de l'objet
// la surface n'a pas de bord (topologique)
template<typename MESH>
void surfaceFacePointMask(const MESH& m,Dart d, std::vector<Dart> &p_point, std::vector<Dart> &q_point)
{
    // parcourt de la face
    Dart t = d;
    do {
        // p-points : sommets de la face
        p_point.push_back(t);
        // q-points : sommets des aretes opposées des faces adjacentes pour chaque coté de la face
        Dart e = surfaceAdjacentFace(m,t);
        q_point.push_back(phi<11>(m,e));
        q_point.push_back(phi_1(m,e));
        // coté suivant
        t = phi1(m,t);
    } while (t != d);
}

// remplit les vecteurs de p/q/r-points pour l'arete incidente à d (d à la surface de l'objet)
// la surface n'a pas de bord (topologique)
template<typename MESH>
void surfaceEdgePointMask(const MESH& m,Dart d,
                          std::vector<Dart> &p_point,
                          std::vector<Dart> &q_point,
                          std::vector<Dart> &r_point)
{
    // extremités de l'arete
    p_point.push_back(d);
    p_point.push_back(phi1(m,d));

    // extremités de l'arete opposée sur la face (en haut)
    q_point.push_back(phi<11>(m,d));
    q_point.push_back(phi_1(m,d));
    // en haut à droite
    Dart t = surfaceAdjacentFace(m,phi1(m,d));
    r_point.push_back(phi_1(m,t));
    // en haut à gauche
    t = surfaceAdjacentFace(m,phi_1(m,d));
    r_point.push_back(phi<11>(m,t));

    // extremités de l'arete opposée sur la face adjacente (en bas)
    Dart e = surfaceAdjacentFace(m,d);
    q_point.push_back(phi<11>(m,e));
    q_point.push_back(phi_1(m,e));
    // en bas à gauche
    t = surfaceAdjacentFace(m,phi1(m,e));
    r_point.push_back(phi_1(m,t));
    // en bas à droite
    t = surfaceAdjacentFace(m,phi_1(m,e));
    r_point.push_back(phi<11>(m,t));
}

/* ----------------------------- BUTTERFLY RULES ----------------------------- */

template<typename MESH,typename T>
Vec3 sum(const MESH& m,const std::vector<Dart> &v,const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>>& attribute)
{
    Vec3 s = Vec3(0.f, 0.f, 0.f);
    for (Dart i : v) s += value<T>(m, attribute,typename MESH::Vertex(i));
    return s;
}

// Le nombre de p-points est nécessairement 8,
// par contre pour les volumes en bordure de l'objet il peut manquer des q-points
// Si c'est le cas, on fait une subdivision linéaire au lieu de la subdivision interpolante
template<typename T,typename MESH>
Vec3 volumePointRule(const MESH& m,
                     const std::vector<Dart> &p_point,
                     const std::vector<Dart> &q_point,
                     const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>>& attribute)
{
    static const float _W_ = 1.f/16;
    if (q_point.size() == 0) {
        return (1.f/8)*sum<MESH>(m,p_point,attribute);
    }
    else {
        return ((6*_W_ + 1)/8)*sum<MESH>(m,p_point,attribute) - (_W_/4)*sum<MESH>(m,q_point,attribute);
    }
}

// Meme principe
template<typename T,typename MESH>
Vec3 facePointRule(const MESH& m,const std::vector<Dart> &p_point, const std::vector<Dart> &q_point,
                   const std::vector<Dart> &r_point, const std::vector<Dart> &s_point,
                   const std::vector<Dart> &t_point,const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>>& attribute)
{
    static const float _W_ = 1.f/16;
    if (q_point.size() == 0) {
        return (1.f/4)*sum<MESH>(m,p_point,attribute);
    }
    else {
        float w1 = (2*_W_ + 1)/4, w2 = _W_/4, w3 = _W_/8;
        return w1*sum<MESH>(m,p_point,attribute) +
                w2*sum<MESH>(m,q_point,attribute) -
                w2*sum<MESH>(m,r_point,attribute) -
                w3*sum<MESH>(m,s_point,attribute) -
                w3*sum<MESH>(m,t_point,attribute);
    }
}

// Même principe, et on prend en plus en considération la possibilité d'avoir
// un nombre de q-points différent de 4*2 (wq toujours > 2)
template<typename T,typename MESH>
Vec3 edgePointRule(const MESH& m,const std::vector<Dart> &p_point, const std::vector<Dart> &q_point,
                   const std::vector<Dart> &r_point, const std::vector<Dart> &s_point,
                   const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>>& attribute)
{
    static const float _W_ = 1.f/16;
    int N = q_point.size()/2;
    float w1 = 1.f/2;
    if (N == 0) {
        return w1*sum<MESH>(m,p_point,attribute);
    }
    else {
        float w2 = _W_/N, w3 = _W_/(2*N);
        return w1*sum<MESH>(m,p_point,attribute) +
                w2*sum<MESH>(m,q_point,attribute) -
                w2*sum<MESH>(m,r_point,attribute) -
                w3*sum<MESH>(m,s_point,attribute);
    }
}

// Le nombre de p/q-points est constant car on ne considère pas les surface à bord
template<typename T,typename MESH>
Vec3 surfaceFacePointRule(const MESH& m,const std::vector<Dart> &p_point,
                          const std::vector<Dart> &q_point,
                          const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>>& attribute)
{
    static const float _W_ = 1.f/16;
    Vec3 sum_p = sum<MESH>(m,p_point,attribute), sum_q = sum<MESH>(m,q_point,attribute);
    float wp = 1.f/4.f + _W_, wq = -_W_/2.f;
    return wp*sum_p + wq*sum_q;
}

// Le nombre de p/q/r-points est constant car on ne considère pas les surface à bord
template<typename T,typename MESH>
Vec3 surfaceEdgePointRule(const MESH& m,const std::vector<Dart> &p_point, const std::vector<Dart> &q_point,
                          const std::vector<Dart> &r_point,const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<T>>& attribute)
{
    static const float _W_ = 1.f/16;
    Vec3 sum_p = sum<MESH>(m,p_point,attribute),
            sum_q = sum<MESH>(m,q_point,attribute),
            sum_r = sum<MESH>(m,r_point,attribute);
    float wp = 1.f/2.f, wq = _W_/2.f, wr = -_W_/2.f;
    return wp*sum_p + wq*sum_q + wr*sum_r;
}


template<typename MESH>
void subdivideEdge(MESH& m,Dart d, Vec3& p,const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Vec3>>& attribute)
{
	cut_edge(m,typename MESH::Edge(d));
	Dart tmp = phi1(m,Dart(0));
	while(Dart(0) != tmp){
		tmp = phi1(m,tmp);
	}
	value<Vec3>(m, attribute,typename MESH::Vertex(phi1(m,d))) = p;
}

template<typename MESH>
Dart subdivideFace(MESH& m,Dart d, Vec3& p,const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Vec3>>& attribute)
{
	using Vertex = typename MESH::Vertex;
	Dart e = phi1(m,d);
	Dart f = phi<1111>(m,e);
	// 1ère coupe
	cut_face(m,Vertex(e),Vertex(f));
	// Nouveau sommet
	subdivideEdge(m,phi_1(m,e),p,attribute);

	Dart result = phi_1(m,e);

	// 2ème et 3ème coupes
	cut_face(m,Vertex(phi_1(m,e)), Vertex(phi<11>(m,e)));
	cut_face(m,Vertex(phi_1(m,f)), Vertex(phi<11>(m,f)));

	return result;
}

// /!\ les faces du volume sont déjà subdivisées
template<typename MESH>
void subdivideVolume(MESH& m,Dart d, Vec3& p,const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Vec3>>& attribute)
{
	using Vertex = typename MESH::Vertex;
	// direction pour la coupe centrale en deux carrés
	Dart first_cut_dir = phi1(m,d);
	// directions pour les coupes des carrés gauches et droits en 2 briques
	Dart left_cut_dir = phi<1>(m,first_cut_dir), right_cut_dir = phi<21112>(m,first_cut_dir);

	Dart cut_direction[4];
	// directions pour les coupes des deux briques du carré gauche
	cut_direction[0] = phi<1211>(m,left_cut_dir);
	cut_direction[1] = phi<1212111>(m,left_cut_dir);
	// directions pour les coupes des deux briques du carré droit
	cut_direction[2] = phi<1211>(m,right_cut_dir);
	cut_direction[3] = phi<1212111>(m,right_cut_dir);
	
	auto cutVolume = [&](Dart d)->typename MESH::Face{
		std::vector<Dart> cut_path;
        Dart t = d;
        do {
            cut_path.push_back(t);
            t = phi<121>(m,t);
        } while (t != d);
        return cut_volume(m,cut_path);
	};

	// première coupe en 2 sous-volumes carrés et subdivision de la face centrale
	typename MESH::Face F = cutVolume(first_cut_dir);
	subdivideFace(m,F.dart, p,attribute);

	// coupes des 2 carrés et coupe de la face entre les 2 briques résultantes pour chaque carré
	F = cutVolume(left_cut_dir);
	cut_face(m,Vertex(F.dart), Vertex(phi<111>(m,F.dart)));
	F = cutVolume(right_cut_dir);
	cut_face(m,Vertex(F.dart), Vertex(phi<111>(m,F.dart)));
	// coupes des briques
	for (int i = 0; i < 4; ++i) {
		cutVolume(cut_direction[i]);
	}
}

template<typename MESH>
void subdivideListEdges(MESH& m,std::vector<Dart>& edges,std::queue<Vec3>& edge_points,const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Vec3>>& attribute){
	for (Dart d : edges) {
		subdivideEdge(m,d, edge_points.front(),attribute);
		edge_points.pop();
	}
}

void subdivideListEdges(MRCmap3& m,std::vector<Dart>& edges,std::queue<Vec3>& edge_points,const std::shared_ptr<typename mesh_traits<MRCmap3>::template Attribute<Vec3>>& attribute){
	for (Dart d : edges) {
		MRCmap3 m2(m,edge_level(m,d)+1);
		subdivideEdge(m2,d, edge_points.front(),attribute);
		edge_points.pop();
	}
}

template<typename MESH>
void subdivideListFaces(MESH& m,std::vector<Dart>& faces,std::queue<Vec3>& face_points,const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Vec3>>& attribute){
	for (Dart d : faces) {
		subdivideFace(m,d, face_points.front(),attribute);
		face_points.pop();
	}
}

void subdivideListFaces(MRCmap3& m,std::vector<Dart>& faces,std::queue<Vec3>& face_points,const std::shared_ptr<typename mesh_traits<MRCmap3>::template Attribute<Vec3>>& attribute){
	for (Dart d : faces) {
		MRCmap3 m2(m,face_level(m,d)+1);
		subdivideFace(m2,d, face_points.front(),attribute);
		face_points.pop();
	}
}

template<typename MESH>
void subdivideListVolumes(MESH& m,std::vector<Dart>& volumes,std::queue<Vec3>& volume_points,const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Vec3>>& attribute){
	for (Dart d : volumes) {
		subdivideVolume(m,d, volume_points.front(),attribute);
		volume_points.pop();
	}
}

void subdivideListVolumes(MRCmap3& m,std::vector<Dart>& volumes,std::queue<Vec3>& volume_points,const std::shared_ptr<typename mesh_traits<MRCmap3>::template Attribute<Vec3>>& attribute){
	for (Dart d : volumes) {
		MRCmap3 m2(m,volume_level(m,d)+1);
		subdivideVolume(m2,d, volume_points.front(),attribute);
		volume_points.pop();
	}
}


template<typename MESH>
void butterflySubdivisionVolumeAdaptative(MESH& m,double angle_threshold,const std::shared_ptr<typename mesh_traits<MESH>::template Attribute<Vec3>>& attribute)
{
	using Volume = typename MESH::Volume;
	using Face = typename MESH::Face;
	using Edge = typename MESH::Edge;
	
	CellMarker<MESH,Edge> cm_surface(m);
	CellMarker<MESH,Volume> cm_cell(m);

	// Sélection de volumes à subdiviser en fonction des angles entre les faces à la surface
	for (Dart d = begin(m); d != end(m);d = next(m,d)) {
		if (is_boundary(m,phi3(m,d)) && !cm_surface.is_marked(Edge(d))){
			// recherche de la seconde face incidente à l'arete en surface
			Dart d2 = phi2(m,d);
			while (!is_boundary(m,phi3(m,d2))) {
				d2 = phi<32>(m,d2);
			}
			auto edge_angle = geometry::angle(m, typename MESH::Face2(d), typename MESH::Face2(d2), attribute.get());
			// Sélection en fonction de l'angle entre les faces
			if (std::abs(edge_angle) > angle_threshold) {
				if (!cm_cell.is_marked(Volume(d)))
					cm_cell.mark(Volume(d));
				if (!cm_cell.is_marked(Volume(d2)))
					cm_cell.mark(Volume(d2));
			}
			cm_surface.mark(Edge(d));
		}
	}

	CellMarker<MESH,Edge>  cm_edge(m);
	CellMarker<MESH,Face>  cm_face(m);
	CellMarker<MESH,Volume>  cm_volume(m);
	std::queue<Vec3> volume_points, face_points, edge_points;
	std::vector<Dart> edges, faces, volumes;
	std::vector<Dart> p_point, q_point, r_point, s_point, t_point;

	// Calcul des plongements pour les points-arete/face/volume
	for (Dart d = begin(m); d!= end(m); d = next(m,d)) {
		if (!is_boundary(m,d) && cm_cell.is_marked(Volume(d))) {
			// pour chaque volume sélectionné, s'il est valide pour la subdivision,
			// alors pour chaque brin du volume, on calcule les points-arete/face/volume
			foreach_dart_of_orbit(m,Volume(d), [&] (Dart t) -> bool{
				// points arete
				if (!cm_edge.is_marked(Edge(t))) {
					if (is_incident_to_boundary(m,Edge(t)) && !is_boundary(m,phi3(m,t))) {
						// on ignore les brins en surface mais en jonction de deux volumes
					}
					// cas surface
					else if (is_boundary(m,phi3(m,t))) {
						p_point.clear();
						q_point.clear();
						r_point.clear();
						surfaceEdgePointMask(m,t, p_point, q_point, r_point);
						
						Vec3 vec = surfaceEdgePointRule<Vec3>(m,p_point, q_point, r_point,attribute);

						edge_points.push(vec);
						edges.push_back(t);
						cm_edge.mark(Edge(t));
					}
					// cas volume
					else {
						p_point.clear();
						q_point.clear();
						r_point.clear();
						s_point.clear();
						edgePointMask(m,t, p_point, q_point, r_point, s_point);
						
						Vec3 vec = edgePointRule<Vec3>(m,p_point, q_point, r_point,s_point,attribute);

						edge_points.push(vec);

						edges.push_back(t);
						cm_edge.mark(Edge(t));
					}
				}
				// points face
				if (!cm_face.is_marked(Face(t))) {
					// cas surface
					if (is_incident_to_boundary(m,Face(t))) {
						p_point.clear();
						q_point.clear();
						surfaceFacePointMask(m,t, p_point, q_point);
						
						Vec3 vec = surfaceFacePointRule<Vec3>(m,p_point, q_point,attribute);

						face_points.push(vec);
					}
					// cas volume
					else {
						p_point.clear();
						q_point.clear();
						r_point.clear();
						s_point.clear();
						t_point.clear();
						facePointMask(m,t, p_point, q_point, r_point, s_point, t_point);

						Vec3 vec = facePointRule<Vec3>(m,p_point, q_point, r_point,s_point,t_point,attribute);
						
						face_points.push(vec);
					}
					faces.push_back(t);
					cm_face.mark(Face(t));
				}
				// points volumes
				if (!cm_volume.is_marked(Volume(t))) {
					p_point.clear();
					q_point.clear();
					volumePointMask(m,t, p_point, q_point);
					
					Vec3 vec = volumePointRule<Vec3>(m,p_point, q_point,attribute);

					volume_points.push(vec);
					volumes.push_back(t);
					cm_volume.mark(Volume(t));
				}
				return true;
			});
			cm_cell.unmark(Volume(d));
		}
	}
	
	MESH m2(m,m.cph().maximum_level()+1);
	
	Dart tes_dart = phi2(m2,Dart(0));
	Dart tmp = phi1(m2,edges[1]);
	while(tmp != edges[1]){
		tmp = phi1(m2,tmp);
	}
	// Subdivision des aretes
	subdivideListEdges(m,edges,edge_points,attribute);
	tmp = phi1(m2,edges[1]);
	while(tmp != edges[1]){
		tmp = phi1(m2,tmp);
	}
	// Subdivision des faces
	subdivideListFaces(m,faces,face_points,attribute);
	// Subdivision des volumes
	subdivideListVolumes(m,volumes,volume_points,attribute);
}

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_SUBDIVISION_H_

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

#ifndef CGOGN_GEOMETRY_ALGOS_EAR_TRIANGULATION_H_
#define CGOGN_GEOMETRY_ALGOS_EAR_TRIANGULATION_H_

#include <cgogn/core/functions/mesh_ops/face.h>

//#include <cgogn/geometry/types/geometry_traits.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/functions/inclusion.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <set>
#define _USE_MATH_DEFINES
#include <math.h>

namespace cgogn
{

namespace geometry
{

template <typename MESH>
class EarTriangulation
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	class VertexPoly
	{
	public:
		using VPMS = std::multiset<VertexPoly*, bool (*)(VertexPoly*, VertexPoly*)>;

		int32 id;
		Vertex vert_;
		Scalar value_;
		Scalar length_;
		VertexPoly* prev_;
		VertexPoly* next_;
		typename VPMS::iterator ear_;

		VertexPoly(Vertex ve, Scalar va, Scalar l, VertexPoly* p = nullptr)
			: vert_(ve), value_(va), length_(l), prev_(p), next_(nullptr)
		{
			if (prev_ != nullptr)
				prev_->next_ = this;
		}

		static void close(VertexPoly* first, VertexPoly* last)
		{
			last->next_ = first;
			first->prev_ = last;
		}

		static VertexPoly* erase(VertexPoly* vp)
		{
			VertexPoly* tmp = vp->prev_;
			tmp->next_ = vp->next_;
			vp->next_->prev_ = tmp;
			delete vp;
			return tmp;
		}
	};

	using VPMS = typename VertexPoly::VPMS;

	// normal to polygon (for orientation of angles)
	Vec3 normalPoly_;

	// ref on map
	MESH& m_;

	// pointer to position attribute
	const typename mesh_traits<MESH>::template Attribute<Vec3>* position_;

	inline const Vec3& POSITION(Vertex v)
	{
		return value<Vec3>(m_, position_, v);
	}

	// map of ears
	VPMS ears_;

	// is current polygon convex
	bool convex_;

	// number of vertices
	uint32 nb_verts_;

	// initial face
	Face face_;

	inline static bool cmp_VP(VertexPoly* lhs, VertexPoly* rhs)
	{
		if (std::abs(lhs->value_ - rhs->value_) < Scalar(0.2))
			return lhs->length_ < rhs->length_;
		return lhs->value_ < rhs->value_;
	}

	void recompute_2_ears(VertexPoly* vp)
	{
		VertexPoly* vprev = vp->prev_;
		VertexPoly* vp2 = vp->next_;
		VertexPoly* vnext = vp2->next_;
		const Vec3& Ta = POSITION(vp->vert_);
		const Vec3& Tb = POSITION(vp2->vert_);
		const Vec3& Tc = POSITION(vprev->vert_);
		const Vec3& Td = POSITION(vnext->vert_);

		// compute angle
		Vec3 v1 = Tb - Ta;
		Vec3 v2 = Tc - Ta;
		Vec3 v3 = Td - Tb;

		v1.normalize();
		v2.normalize();
		v3.normalize();

		Scalar dotpr1 = std::acos(v1.dot(v2)) / Scalar(M_PI_2);
		Scalar dotpr2 = std::acos(-(v1.dot(v3))) / Scalar(M_PI_2);

		if (!convex_) // if convex no need to test if vertex is an ear (yes)
		{
			Vec3 nv1 = v1.cross(v2);
			Vec3 nv2 = v1.cross(v3);

			if (nv1.dot(normalPoly_) < Scalar(0))
				dotpr1 = Scalar(10) - dotpr1; // not an ear (concave)
			if (nv2.dot(normalPoly_) < Scalar(0))
				dotpr2 = Scalar(10) - dotpr2; // not an ear (concave)

			bool finished = (dotpr1 >= Scalar(5)) && (dotpr2 >= Scalar(5));
			for (auto it = ears_.rbegin(); (!finished) && (it != ears_.rend()) && ((*it)->value_ > Scalar(5)); ++it)
			{
				Vertex id = (*it)->vert_;
				const Vec3& P = POSITION(id);

				if ((dotpr1 < Scalar(5)) && (id.dart != vprev->vert_.dart))
					if (in_triangle(P, normalPoly_, Tb, Tc, Ta))
						dotpr1 = Scalar(5); // not an ear !

				if ((dotpr2 < Scalar(5)) && (id.dart != vnext->vert_.dart))
					if (in_triangle(P, normalPoly_, Td, Ta, Tb))
						dotpr2 = Scalar(5); // not an ear !

				finished = (dotpr1 >= Scalar(5)) && (dotpr2 >= Scalar(5));
			}
		}

		vp->value_ = dotpr1;
		vp->length_ = Scalar((Tb - Tc).squaredNorm());
		vp->ear_ = ears_.insert(vp);
		vp2->value_ = dotpr2;
		vp->length_ = Scalar((Td - Ta).squaredNorm());
		vp2->ear_ = ears_.insert(vp2);

		// polygon if convex only if all vertices have convex angle (last have ...)
		convex_ = (*(ears_.rbegin()))->value_ < Scalar(5);
	}

	Scalar ear_angle(const Vec3& P1, const Vec3& P2, const Vec3& P3)
	{
		Vec3 v1 = P1 - P2;
		Vec3 v2 = P3 - P2;
		v1.normalize();
		v2.normalize();

		Scalar dotpr = std::acos(v1.dot(v2)) / Scalar(M_PI_2);

		Vec3 vn = v1.cross(v2);
		if (vn.dot(normalPoly_) > Scalar(0))
			dotpr = Scalar(10) - dotpr; // not an ear (concave, store at the end for optimized use for intersections)

		return dotpr;
	}

	bool ear_intersection(VertexPoly* vp)
	{
		VertexPoly* endV = vp->prev_;
		VertexPoly* curr = vp->next_;
		const Vec3& Ta = POSITION(vp->vert_);
		const Vec3& Tb = POSITION(curr->vert_);
		const Vec3& Tc = POSITION(endV->vert_);
		curr = curr->next_;

		while (curr != endV)
		{
			if (in_triangle(POSITION(curr->vert_), normalPoly_, Tb, Tc, Ta))
			{
				vp->value_ = Scalar(5); // not an ear !
				return false;
			}
			curr = curr->next_;
		}
		return true;
	}

	////////////////////////////////
	// CMapBase (and convertible) //
	////////////////////////////////

	template <typename MESHTYPE, typename std::enable_if_t<std::is_convertible_v<MESHTYPE&, CMapBase&>>* = nullptr>
	std::tuple<VertexPoly*, VertexPoly*, uint32, bool> init_chained_vertexpoly_list(
		MESHTYPE& m, const typename mesh_traits<MESHTYPE>::Face f,
		const typename mesh_traits<MESHTYPE>::template Attribute<Vec3>* position)
	{
		VertexPoly* vpp = nullptr;
		VertexPoly* prem = nullptr;
		uint32 nb_verts = 0;
		bool convex = true;

		Dart a = f.dart;
		Dart b = phi1(m_, a);
		Dart c = phi1(m_, b);
		do
		{
			const Vec3& P1 = POSITION(Vertex(a));
			const Vec3& P2 = POSITION(Vertex(b));
			const Vec3& P3 = POSITION(Vertex(c));

			Scalar val = ear_angle(P1, P2, P3);
			VertexPoly* vp = new VertexPoly(Vertex(b), val, Scalar((P3 - P1).squaredNorm()), vpp);

			if (vp->value_ > Scalar(5)) // concav angle
				convex = false;

			if (vpp == nullptr)
				prem = vp;
			vpp = vp;
			a = b;
			b = c;
			c = phi1(m_, c);
			nb_verts++;
		} while (a != f.dart);

		return {vpp, prem, nb_verts, convex};
	}

public:
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(EarTriangulation);

	/**
	 * @brief EarTriangulation constructor
	 * @param map ref on map
	 * @param f the face to tringulate
	 * @param position attribute of position to use
	 */
	EarTriangulation(MESH& mesh, const typename mesh_traits<MESH>::Face f,
					 const typename mesh_traits<MESH>::template Attribute<Vec3>* position)
		: m_(mesh), position_(position), ears_(cmp_VP)
	{
		uint32 codeg = codegree(m_, f);

		if (codeg <= 3)
		{
			face_ = f;
			nb_verts_ = codeg;
			return;
		}

		// compute normals for orientation
		normalPoly_ = normal(m_, f, position_);

		// first pass create polygon in chained list with angle computation
		auto [vpp, prem, nb_verts, convex] = init_chained_vertexpoly_list(mesh, f, position);
		nb_verts_ = nb_verts;
		convex_ = convex;

		VertexPoly::close(prem, vpp);

		if (convex_)
		{
			// second pass with no test of intersections with polygons
			vpp = prem;
			for (uint32 i = 0; i < nb_verts_; ++i)
			{
				vpp->ear_ = ears_.insert(vpp);
				vpp = vpp->next_;
			}
		}
		else
		{
			// second pass test intersections with polygons
			vpp = prem;
			for (uint32 i = 0; i < nb_verts_; ++i)
			{
				if (vpp->value_ < Scalar(5))
					ear_intersection(vpp);
				vpp->ear_ = ears_.insert(vpp);
				vpp = vpp->next_;
			}
		}
	}

	/**
	 * @brief compute table of vertices indices (embeddings) of triangulation
	 * @param table_indices
	 */
	template <typename FUNC>
	void append_indices(std::vector<uint32>& table_indices, const FUNC& post_func)
	{
		if (nb_verts_ < 3)
			return;

		//		table_indices.reserve((nb_verts_ - 2) * 3);

		if (nb_verts_ == 3)
		{
			foreach_incident_vertex(m_, face_, [&](Vertex v) -> bool {
				table_indices.push_back(index_of(m_, v));
				return true;
			});
			post_func();
			return;
		}

		while (nb_verts_ > 3)
		{
			// take best (and valid!) ear
			typename VPMS::iterator be_it = ears_.begin(); // best ear
			VertexPoly* be = *be_it;

			table_indices.push_back(index_of(m_, be->vert_));
			table_indices.push_back(index_of(m_, be->next_->vert_));
			table_indices.push_back(index_of(m_, be->prev_->vert_));
			post_func();
			--nb_verts_;

			if (nb_verts_ > 3) // do not recompute if only one triangle left
			{
				// remove ears and two sided ears
				ears_.erase(be_it); // from map of ears
				ears_.erase(be->next_->ear_);
				ears_.erase(be->prev_->ear_);
				be = VertexPoly::erase(be); // and remove ear vertex from polygon
				recompute_2_ears(be);
			}
			else // finish (no need to update ears)
			{
				// remove ear from polygon
				be = VertexPoly::erase(be);
				// last triangle
				table_indices.push_back(index_of(m_, be->vert_));
				table_indices.push_back(index_of(m_, be->next_->vert_));
				table_indices.push_back(index_of(m_, be->prev_->vert_));
				post_func();
				// release memory of last triangle in polygon
				delete be->next_;
				delete be->prev_;
				delete be;
			}
		}
	}

	/**
	 * @brief apply the ear triangulation the face
	 */
	void apply()
	{
		while (nb_verts_ > 3)
		{
			// take best (and valid!) ear
			typename VPMS::iterator be_it = ears_.begin(); // best ear
			VertexPoly* be = *be_it;

			cut_face(m_, be->prev_->vert_, be->next_->vert_);

			--nb_verts_;

			if (nb_verts_ > 3) // do not recompute if only one triangle left
			{
				// remove ears and two sided ears
				ears_.erase(be_it); // from map of ears
				ears_.erase(be->next_->ear_);
				ears_.erase(be->prev_->ear_);
				// replace dart to be in remaining poly
				be->prev_->vert_ = Vertex(phi2(m_, phi_1(m_, be->prev_->vert_.dart)));
				be = VertexPoly::erase(be); // and remove ear vertex from polygon
				recompute_2_ears(be);
			}
			else // finish (no need to update ears)
			{
				// release memory of last triangle in polygon
				delete be->next_;
				delete be->prev_;
				delete be;
			}
		}
	}
};

/**
 * @brief compute ear triangulation
 * @param map
 * @param f face
 * @param position
 * @param table_indices table of indices (vertex embedding) to append
 */
template <typename MESH, typename FUNC>
void append_ear_triangulation(const MESH& mesh, const typename mesh_traits<MESH>::Face f,
							  const typename mesh_traits<MESH>::template Attribute<Vec3>* position,
							  std::vector<uint32>& table_indices, const FUNC& post_func)
{
	EarTriangulation tri(const_cast<MESH&>(mesh), f, position);
	tri.append_indices(table_indices, post_func);
}

/**
 * @brief apply ear triangulation to a face (face is cut)
 * @param map
 * @param f
 * @param position
 */
template <typename MESH>
void apply_ear_triangulation(MESH& mesh, const typename mesh_traits<MESH>::Face f,
							 const typename mesh_traits<MESH>::template Attribute<Vec3>* position)
{
	EarTriangulation tri(mesh, f, position);
	tri.apply();
}

/**
 * @brief apply ear triangulation to a map
 * @param map
 * @param f
 * @param position
 */
template <typename MESH>
void apply_ear_triangulation(MESH& mesh, const typename mesh_traits<MESH>::template Attribute<Vec3>* position)
{
	foreach_cell(mesh, [&](typename mesh_traits<MESH>::Face f) -> bool {
		if (codegree(mesh, f) > 3)
		{
			EarTriangulation tri(mesh, f, position);
			tri.apply();
		}
		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_EAR_TRIANGULATION_H_

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

#ifndef CGOGN_MODELING_DECIMATION_EDGE_TRAVERSOR_EDGE_LENGTH_H_
#define CGOGN_MODELING_DECIMATION_EDGE_TRAVERSOR_EDGE_LENGTH_H_

#include <cgogn/modeling/algos/decimation/edge_traversor.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/attributes.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/algos/length.h>

#include <map>

namespace cgogn
{

namespace modeling
{

////////////////////////////////////////
// EdgeTraversor_EdgeLength CRTP Base //
////////////////////////////////////////

template <typename Derived, typename MESH, typename VEC>
class EdgeTraversor_EdgeLength_Base : public EdgeTraversor<MESH>
{
public:

	using Self = EdgeTraversor_EdgeLength_Base<Derived, MESH, VEC>;

	using Scalar = typename geometry::vector_traits<VEC>::Scalar;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	struct EdgeInfo
	{
		typename std::multimap<Scalar, Edge>::const_iterator it_;
		bool valid_;

		EdgeInfo() : valid_(false) {}
	};

protected:

	const typename mesh_traits<MESH>::template AttributePtr<VEC> vertex_position_;
	typename mesh_traits<MESH>::template AttributePtr<EdgeInfo> einfo_;
	std::multimap<Scalar, Edge> edges_;
	Edge e1_, e2_;

public:

	EdgeTraversor_EdgeLength_Base(
		MESH& m,
		const typename mesh_traits<MESH>::template AttributePtr<VEC> vertex_position
	) : EdgeTraversor<MESH>(m),
		vertex_position_(vertex_position)
	{
		einfo_ = add_attribute<EdgeInfo, Edge>(this->m_, "EdgeTraversor_EdgeLength_EdgeInfo");

		foreach_cell(this->m_, [&] (Edge e) -> bool
		{
			value<EdgeInfo>(this->m_, einfo_, e).valid_ = false;
			update_edge_info(e);
			return true;
		});
	}

	~EdgeTraversor_EdgeLength_Base() override
	{
		remove_attribute<Edge>(this->m_, einfo_);
	}

	void update_edge_info(Edge e)
	{
		EdgeInfo& ei = value<EdgeInfo>(this->m_, einfo_, e);
		if (edge_can_collapse(this->m_, e))
		{
			Scalar cost = geometry::length<VEC>(this->m_, e, vertex_position_);
			if (ei.valid_)
			{
				edges_.erase(ei.it_);
				ei.it_ = edges_.insert(std::make_pair(cost, e));
			}
			else
			{
				ei.it_ = edges_.insert(std::make_pair(cost, e));
				ei.valid_ = true;
			}
		}
		else
		{
			if (ei.valid_)
			{
				edges_.erase(ei.it_);
				ei.valid_ = false;
			}
		}
	}

	class const_iterator
	{
	public:

		const Self* const trav_ptr_;
		typename std::multimap<Scalar, Edge>::const_iterator edge_it_;

		inline const_iterator(const EdgeTraversor_EdgeLength_Base<Derived, MESH, VEC>* trav, typename std::multimap<Scalar, Edge>::const_iterator it) :
			trav_ptr_(trav),
			edge_it_(it)
		{}

		inline const_iterator(const const_iterator& it) :
			trav_ptr_(it.trav_ptr_),
			edge_it_(it.edge_it_)
		{}

		inline const_iterator& operator=(const const_iterator& it)
		{
			trav_ptr_ = it.trav_ptr_;
			edge_it_ = it.edge_it_;
			return *this;
		}

		inline const_iterator& operator++()
		{
			edge_it_ = trav_ptr_->edges_.begin();
			return *this;
		}

		inline const Edge& operator*() const
		{
			return (*edge_it_).second;
		}

		inline bool operator!=(const_iterator it) const
		{
			cgogn_assert(trav_ptr_ == it.trav_ptr_);
			return edge_it_ != it.edge_it_;
		}
	};

	const_iterator begin() const
	{
		return const_iterator(this, edges_.begin());
	}

	const_iterator end() const
	{
		return const_iterator(this, edges_.end());
	}
};

//////////////////////////////
// EdgeTraversor_EdgeLength //
//////////////////////////////

template <typename MESH, typename VEC>
class EdgeTraversor_EdgeLength : public EdgeTraversor_EdgeLength_Base<EdgeTraversor_EdgeLength<MESH, VEC>, MESH, VEC>
{

};

////////////////////////////////////
// EdgeTraversor_EdgeLength CMap2 //
////////////////////////////////////

template <typename VEC>
class EdgeTraversor_EdgeLength<CMap2, VEC> : public EdgeTraversor_EdgeLength_Base<EdgeTraversor_EdgeLength<CMap2, VEC>, CMap2, VEC>
{
public:

	using Inherit = EdgeTraversor_EdgeLength_Base<EdgeTraversor_EdgeLength<CMap2, VEC>, CMap2, VEC>;
	using EdgeInfo = typename Inherit::EdgeInfo;

	using Scalar = typename geometry::vector_traits<VEC>::Scalar;
	using Vertex = typename mesh_traits<CMap2>::Vertex;
	using Edge = typename mesh_traits<CMap2>::Edge;

	EdgeTraversor_EdgeLength(
		CMap2& m,
		const typename mesh_traits<CMap2>::template AttributePtr<VEC> vertex_position
	) : Inherit(m, vertex_position)
	{}

	~EdgeTraversor_EdgeLength() override
	{}

	void pre_collapse(Edge e) override
	{
		EdgeInfo& ei = value<EdgeInfo>(this->m_, this->einfo_, e);
		if (ei.valid_) { this->edges_.erase(ei.it_); ei.valid_ = false; }

		this->e1_ = Edge(this->m_.phi2(this->m_.phi_1(e.dart)));
		this->e2_ = Edge(this->m_.phi2(this->m_.phi_1(this->m_.phi2(e.dart))));

		Dart ed1 = e.dart;
		Dart ed2 = this->m_.phi2(ed1);

		EdgeInfo& ei1 = value<EdgeInfo>(this->m_, this->einfo_, Edge(this->m_.phi1(ed1)));
		if (ei1.valid_) { this->edges_.erase(ei1.it_); ei1.valid_ = false; }

		EdgeInfo& ei2 = value<EdgeInfo>(this->m_, this->einfo_, Edge(this->m_.phi_1(ed1)));
		if (ei2.valid_) { this->edges_.erase(ei2.it_); ei2.valid_ = false; }

		EdgeInfo& ei3 = value<EdgeInfo>(this->m_, this->einfo_, Edge(this->m_.phi1(ed2)));
		if (ei3.valid_) { this->edges_.erase(ei3.it_); ei3.valid_ = false; }

		EdgeInfo& ei4 = value<EdgeInfo>(this->m_, this->einfo_, Edge(this->m_.phi_1(ed2)));
		if (ei4.valid_) { this->edges_.erase(ei4.it_); ei4.valid_ = false; }
	}

	void post_collapse() override
	{
		Dart vit = this->e1_.dart;
		do
		{
			this->update_edge_info(Edge(this->m_.phi1(vit)));
			if(vit == this->e1_.dart || vit == this->e2_.dart)
			{
				this->update_edge_info(Edge(vit));

				Dart vit2 = this->m_.template phi<121>(vit);
				Dart stop = this->m_.phi2(vit);
				do
				{
					this->update_edge_info(Edge(vit2));
					this->update_edge_info(Edge(this->m_.phi1(vit2)));
					vit2 = this->m_.phi1(this->m_.phi2(vit2));
				} while(vit2 != stop);
			}
			else
				this->update_edge_info(Edge(vit));

			vit = this->m_.phi2(this->m_.phi_1(vit));
		} while(vit != this->e1_.dart);
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_DECIMATION_EDGE_TRAVERSOR_EDGE_LENGTH_H_

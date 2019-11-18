/*******************************************************************************
* CGoGN                                                                        *
* Copyright (C) 2019, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_MODULE_SHALLOW_WATER_H_
#define CGOGN_MODULE_SHALLOW_WATER_H_

#include <cgogn/ui/module.h>

#include <cgogn/core/types/cell_marker.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/halfedge.h>
#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/length.h>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class ShallowWater : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension == 2, "ShallowWater can only be used with meshes of dimension 2");

    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

    using Vertex = typename mesh_traits<MESH>::Vertex;
    using HalfEdge = typename mesh_traits<MESH>::HalfEdge;
    using Edge = typename mesh_traits<MESH>::Edge;
    using Face = typename mesh_traits<MESH>::Face;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

	enum BoundaryCondition
	{
		WALL = 0
	};

public:

	ShallowWater(const App& app) :
		ViewModule(app, "ShallowWater (" + std::string{mesh_traits<MESH>::name} + ")")
	{}
	~ShallowWater()
	{
		if (edge_left_side_)
			delete edge_left_side_;
	}

protected:

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));

		domain2D_ = mesh_provider_->load_surface_from_file("/Users/kraemer/Data/surface/grid_tri.off");

		edge_left_side_ = new CellMarker<MESH, HalfEdge>(*domain2D_);

		vertex_position_ = get_attribute<Vec3, Vertex>(*domain2D_, "position");
		vertex_water_position_ = add_attribute<Vec3, Vertex>(*domain2D_, "water_position");

		face_phi_ = add_attribute<Scalar, Face>(*domain2D_, "phi");
		face_zb_ = add_attribute<Scalar, Face>(*domain2D_, "zb");
		face_h_ = add_attribute<Scalar, Face>(*domain2D_, "h");
		face_q_ = add_attribute<Scalar, Face>(*domain2D_, "q");
		face_r_ = add_attribute<Scalar, Face>(*domain2D_, "r");
		face_centroid_ = add_attribute<Vec3, Face>(*domain2D_, "centroid");
		face_area_ = add_attribute<Scalar, Face>(*domain2D_, "area");
		face_swept_ = add_attribute<Scalar, Face>(*domain2D_, "swept");
		face_discharge_ = add_attribute<Scalar, Face>(*domain2D_, "discharge");

		edge_f1_ = add_attribute<Scalar, Edge>(*domain2D_, "f1");
		edge_f2_ = add_attribute<Scalar, Edge>(*domain2D_, "f2");
		edge_f3_ = add_attribute<Scalar, Edge>(*domain2D_, "f3");
		edge_s2L_ = add_attribute<Scalar, Edge>(*domain2D_, "s2L");
		edge_s2R_ = add_attribute<Scalar, Edge>(*domain2D_, "s2R");
		edge_normX_ = add_attribute<Scalar, Edge>(*domain2D_, "normX");
		edge_normY_ = add_attribute<Scalar, Edge>(*domain2D_, "normY");
		edge_length_ = add_attribute<Scalar, Edge>(*domain2D_, "length");
		edge_bc_value_ = add_attribute<Scalar, Edge>(*domain2D_, "bc_value");
		edge_bc_type_ = add_attribute<BoundaryCondition, Edge>(*domain2D_, "bc_type");

		vertex_water_position_->copy(vertex_position_.get());

		parallel_foreach_cell(*domain2D_, [&] (Edge e) -> bool
		{
			if (is_incident_to_boundary(*domain2D_, e))
				value<BoundaryCondition>(*domain2D_, edge_bc_type_, e) = WALL;
			
			Scalar l = geometry::length(*domain2D_, e, vertex_position_.get());
			std::vector<Vertex> vertices = incident_vertices(*domain2D_, e);
			Vec3 vec = value<Vec3>(*domain2D_, vertex_position_, vertices[1]) - value<Vec3>(*domain2D_, vertex_position_, vertices[0]);
			value<Scalar>(*domain2D_, edge_length_, e) = l;
			value<Scalar>(*domain2D_, edge_normX_, e) = vec[1] / l;
			value<Scalar>(*domain2D_, edge_normY_, e) = -vec[0] / l;
			
			std::vector<HalfEdge> halfedges = incident_halfedges(*domain2D_, e);
			edge_left_side_->mark(halfedges[0]);
			edge_left_side_->unmark(halfedges[1]);

			return true;
		});

		parallel_foreach_cell(*domain2D_, [&] (Face f) -> bool
		{
			uint32 nbv = 0;
			Scalar zb = 0.0;
			Vec3 centroid{0, 0, 0};
			foreach_incident_vertex(*domain2D_, f, [&] (Vertex v) -> bool
			{
				const Vec3& p = value<Vec3>(*domain2D_, vertex_position_, v);
				zb += p[2];
				centroid += p;
				nbv++;
				return true;
			});
			zb /= nbv;
			centroid /= nbv;
			Scalar h;
			if (centroid[0] < 0.5) h = 0.5 - centroid[0];
			else h = 0.1;
			
			value<Scalar>(*domain2D_, face_phi_, f) = phi_default_;
			value<Scalar>(*domain2D_, face_zb_, f) = zb;
			value<Scalar>(*domain2D_, face_h_, f) = h;
			value<Scalar>(*domain2D_, face_q_, f) = 0.0;
			value<Scalar>(*domain2D_, face_r_, f) = 0.0;
			value<Vec3>(*domain2D_, face_centroid_, f) = centroid;
			value<Scalar>(*domain2D_, face_area_, f) = geometry::area(*domain2D_, f, vertex_position_.get());

			return true;
		});

		running_ = false;
	}

	void start()
	{
		running_ = true;
	}

	void step()
	{

	}

	void stop()
	{
		running_ = false;
	}

	void update_render_data()
	{

	}

    void interface() override
	{
		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});

		if (!running_)
		{
			if (ImGui::Button("Start"))
				start();
			if (ImGui::Button("step"))
				step();
		}
		else
		{
			if (ImGui::Button("Stop"))
				stop();
		}
		
		
		ImGui::End();
	}

private:

	MeshProvider<MESH>* mesh_provider_;
	MESH* domain2D_;

	CellMarker<MESH, HalfEdge>* edge_left_side_;

	std::shared_ptr<Attribute<Vec3>> vertex_position_;
	std::shared_ptr<Attribute<Vec3>> vertex_water_position_;

	std::shared_ptr<Attribute<Scalar>> face_phi_;
	std::shared_ptr<Attribute<Scalar>> face_zb_;
	std::shared_ptr<Attribute<Scalar>> face_h_;
	std::shared_ptr<Attribute<Scalar>> face_q_;
	std::shared_ptr<Attribute<Scalar>> face_r_;
	std::shared_ptr<Attribute<Vec3>> face_centroid_;
	std::shared_ptr<Attribute<Scalar>> face_area_;
	std::shared_ptr<Attribute<Scalar>> face_swept_;
	std::shared_ptr<Attribute<Scalar>> face_discharge_;

	std::shared_ptr<Attribute<Scalar>> edge_f1_;
	std::shared_ptr<Attribute<Scalar>> edge_f2_;
	std::shared_ptr<Attribute<Scalar>> edge_f3_;
	std::shared_ptr<Attribute<Scalar>> edge_s2L_;
	std::shared_ptr<Attribute<Scalar>> edge_s2R_;
	std::shared_ptr<Attribute<Scalar>> edge_normX_;
	std::shared_ptr<Attribute<Scalar>> edge_normY_;
	std::shared_ptr<Attribute<Scalar>> edge_length_;
	std::shared_ptr<Attribute<Scalar>> edge_bc_value_;
	std::shared_ptr<Attribute<BoundaryCondition>> edge_bc_type_;

	Scalar phi_default_ = 1.0;
	Scalar kx_ = 35.0;
	Scalar ky_ = 35.0;
	Scalar alphaK_ = 0.0;

	bool running_ = false;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SHALLOW_WATER_H_

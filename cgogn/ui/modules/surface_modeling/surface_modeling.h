/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_MODULE_SURFACE_MODELING_H_
#define CGOGN_MODULE_SURFACE_MODELING_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/modeling/algos/decimation/decimation.h>
#include <cgogn/modeling/algos/mesh_repair.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/modeling/algos/topstoc.h>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceModeling : public Module
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "SurfaceModeling can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	SurfaceModeling(const App& app)
		: Module(app, "SurfaceModeling (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  selected_vertex_position_(nullptr)
	{
	}
	~SurfaceModeling()
	{
	}

	void fill_holes(MESH& m)
	{
		modeling::fill_holes(m);
		mesh_provider_->emit_connectivity_changed(&m);
	}

	void remove_small_components(MESH& m, uint32 min_vertices)
	{
		modeling::remove_small_components(m, min_vertices);
		mesh_provider_->emit_connectivity_changed(&m);
	}

	void reverse_orientation(MESH& m)
	{
		cgogn::reverse_orientation(m);
		mesh_provider_->emit_connectivity_changed(&m);
	}

	void triangulate_mesh(MESH& m, Attribute<Vec3>* vertex_position)
	{
		geometry::apply_ear_triangulation(m, vertex_position);
		mesh_provider_->emit_connectivity_changed(&m);
	}

	void decimate_mesh(MESH& m, Attribute<Vec3>* vertex_position)
	{
		modeling::decimate(m, vertex_position, mesh_provider_->mesh_data(&m)->template nb_cells<Vertex>() / 10);
		mesh_provider_->emit_connectivity_changed(&m);
		mesh_provider_->emit_attribute_changed(&m, vertex_position);
	}

	void simplify_mesh(MESH& m, Attribute<Vec3>* vertex_position)
	{
		modeling::topstoc(mesh_provider_, m, vertex_position,
						  0.5 * mesh_provider_->mesh_data(&m)->template nb_cells<Vertex>());
		mesh_provider_->emit_connectivity_changed(&m);
		mesh_provider_->emit_attribute_changed(&m, vertex_position);
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
	}

	void interface() override
	{
		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Surface", [&](MESH* m) {
			selected_mesh_ = m;
			selected_vertex_position_.reset();
			mesh_provider_->mesh_data(selected_mesh_)->outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_position_, "Position",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position_ = attribute; });

			if (selected_vertex_position_)
			{
				if (ImGui::Button("Triangulate"))
					triangulate_mesh(*selected_mesh_, selected_vertex_position_.get());
				if (ImGui::Button("Fill holes"))
					fill_holes(*selected_mesh_);
				static int32 min_vertices = 1000;
				ImGui::SliderInt("Min nb vertices", &min_vertices, 1, 10000);
				if (ImGui::Button("Remove small components"))
					remove_small_components(*selected_mesh_, uint32(min_vertices));
				if (ImGui::Button("Reverse orientation"))
					reverse_orientation(*selected_mesh_);
				if (ImGui::Button("Decimate"))
					decimate_mesh(*selected_mesh_, selected_vertex_position_.get());
				if (ImGui::Button("Simplify"))
					simplify_mesh(*selected_mesh_, selected_vertex_position_.get());
			}
		}
	}

private:
	MESH* selected_mesh_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_MODELING_H_

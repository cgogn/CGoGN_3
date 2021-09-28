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

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/modeling/algos/decimation/decimation.h>
#include <cgogn/modeling/algos/mesh_repair.h>
#include <cgogn/modeling/algos/remeshing/pliant_remeshing.h>
#include <cgogn/modeling/algos/remeshing/topstoc.h>
#include <cgogn/modeling/algos/subdivision.h>

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
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

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
		mesh_provider_->emit_connectivity_changed(m);
	}

	void remove_small_components(MESH& m, uint32 min_vertices)
	{
		modeling::remove_small_components(m, min_vertices);
		mesh_provider_->emit_connectivity_changed(m);
	}

	void reverse_orientation(MESH& m)
	{
		cgogn::reverse_orientation(m);
		mesh_provider_->emit_connectivity_changed(m);
	}

	void triangulate_mesh(MESH& m, Attribute<Vec3>* vertex_position)
	{
		geometry::apply_ear_triangulation(m, vertex_position);
		mesh_provider_->emit_connectivity_changed(m);
	}

	void quadrangulate_mesh(MESH& m, Attribute<Vec3>* vertex_position)
	{
		cgogn::modeling::quadrangulate_all_faces(
			m,
			[&](Vertex v) {
				std::vector<Vertex> av = adjacent_vertices_through_edge(m, v);
				cgogn::value<Vec3>(m, vertex_position, v) = 0.5 * (cgogn::value<Vec3>(m, vertex_position, av[0]) +
																   cgogn::value<Vec3>(m, vertex_position, av[1]));
			},
			[&](Vertex v) {
				Vec3 center;
				center.setZero();
				uint32 count = 0;
				foreach_adjacent_vertex_through_edge(m, v, [&](Vertex av) -> bool {
					center += cgogn::value<Vec3>(m, vertex_position, av);
					++count;
					return true;
				});
				center /= Scalar(count);
				cgogn::value<Vec3>(m, vertex_position, v) = center;
			});

		mesh_provider_->emit_connectivity_changed(m);
		mesh_provider_->emit_attribute_changed(m, vertex_position);
	}

	void delaunay_flips(MESH& m, Attribute<Vec3>* vertex_position)
	{
		foreach_cell(m, [&](Edge e) -> bool {
			std::vector<CMap2::Vertex> iv = incident_vertices(m, e);
			if (degree(m, iv[0]) < 4 || degree(m, iv[1]) < 4)
				return true;
			std::vector<Scalar> op_angles = geometry::opposite_angles(m, e, vertex_position);
			if (op_angles[0] + op_angles[1] > M_PI)
				flip_edge(m, e);
			return true;
		});

		mesh_provider_->emit_connectivity_changed(m);
	}

	void decimate_mesh(MESH& m, Attribute<Vec3>* vertex_position, uint32 percent_vertices_to_remove)
	{
		modeling::decimate(m, vertex_position,
						   mesh_provider_->mesh_data(m).template nb_cells<Vertex>() *
							   (percent_vertices_to_remove / 100.0));
		mesh_provider_->emit_connectivity_changed(m);
		mesh_provider_->emit_attribute_changed(m, vertex_position);
	}

	void simplify_mesh(MESH& m, Attribute<Vec3>* vertex_position)
	{
		modeling::topstoc(mesh_provider_, m, vertex_position,
						  0.5 * mesh_provider_->mesh_data(m).template nb_cells<Vertex>());
		mesh_provider_->emit_connectivity_changed(m);
		mesh_provider_->emit_attribute_changed(m, vertex_position);
	}

	void remesh(MESH& m, Attribute<Vec3>* vertex_position, Scalar edge_length_ratio, bool preserve_features)
	{
		modeling::pliant_remeshing(m, vertex_position, edge_length_ratio, preserve_features);
		mesh_provider_->emit_connectivity_changed(m);
		mesh_provider_->emit_attribute_changed(m, vertex_position);
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
	}

	void interface() override
	{
		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Surface", [&](MESH& m) {
			selected_mesh_ = &m;
			selected_vertex_position_.reset();
			// selected_vertex_normal_.reset();
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_position_, "Position",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position_ = attribute; });

			// imgui_combo_attribute<Vertex, Vec3>(
			// 	*selected_mesh_, selected_vertex_normal_, "Normal",
			// 	[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_normal_ = attribute; });

			if (selected_vertex_position_)
			{
				if (ImGui::Button("Triangulate"))
					triangulate_mesh(*selected_mesh_, selected_vertex_position_.get());
				if (ImGui::Button("Quadrangulate"))
					quadrangulate_mesh(*selected_mesh_, selected_vertex_position_.get());
				if (ImGui::Button("Delaunay flips"))
					delaunay_flips(*selected_mesh_, selected_vertex_position_.get());
				if (ImGui::Button("Fill holes"))
					fill_holes(*selected_mesh_);
				static int32 min_vertices = 1000;
				ImGui::SliderInt("Min nb vertices", &min_vertices, 1, 10000);
				if (ImGui::Button("Remove small components"))
					remove_small_components(*selected_mesh_, uint32(min_vertices));
				if (ImGui::Button("Reverse orientation"))
					reverse_orientation(*selected_mesh_);
				static int32 percent_vertices_to_keep = 90;
				ImGui::SliderInt("% vertices to keep", &percent_vertices_to_keep, 1, 99);
				if (ImGui::Button("Decimate"))
					decimate_mesh(*selected_mesh_, selected_vertex_position_.get(), 100 - percent_vertices_to_keep);
				if (ImGui::Button("Simplify"))
					simplify_mesh(*selected_mesh_, selected_vertex_position_.get());
				static float remesh_edge_length_ratio = 1.0f;
				ImGui::SliderFloat("Edge length target w.r.t. mean", &remesh_edge_length_ratio, 0.0, 3.0);
				static bool preserve_features = false;
				ImGui::Checkbox("Preserve features", &preserve_features);
				if (ImGui::Button("Remesh"))
					remesh(*selected_mesh_, selected_vertex_position_.get(), remesh_edge_length_ratio,
						   preserve_features);
			}
		}
	}

private:
	MESH* selected_mesh_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	// std::shared_ptr<Attribute<Vec3>> selected_vertex_normal_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_MODELING_H_

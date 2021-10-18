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

#ifndef CGOGN_MODULE_SURFACE_HEAT_METHOD_H_
#define CGOGN_MODULE_SURFACE_HEAT_METHOD_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceHeatMethod : public Module
{
	static_assert(mesh_traits<MESH>::dimension == 2, "SurfaceHeatMethod can only be used with meshes of dimension 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	SurfaceHeatMethod(const App& app)
		: Module(app, "SurfaceHeatMethod (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  selected_vertex_position_(nullptr), selected_vertices_set_(nullptr)
	{
	}
	~SurfaceHeatMethod()
	{
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
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);

			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_position_, "Position",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position_ = attribute; });

			imgui_combo_cells_set(md, selected_vertices_set_, "Source vertices",
								  [&](CellsSet<MESH, Vertex>* cs) { selected_vertices_set_ = cs; });

			if (selected_vertex_position_)
			{
			}
		}
	}

private:
	MESH* selected_mesh_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	CellsSet<MESH, Vertex>* selected_vertices_set_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_HEAT_METHOD_H_

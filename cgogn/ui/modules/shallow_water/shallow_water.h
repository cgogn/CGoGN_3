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

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>

#include <cgogn/simulation/algos/shallow_water/shallow_water.h>

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
    using Face = typename mesh_traits<MESH>::Face;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

public:

	ShallowWater(const App& app) :
		ViewModule(app, "ShallowWater (" + std::string{mesh_traits<MESH>::name} + ")")
	{}
	~ShallowWater()
	{}

protected:

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));

		domain2D_ = mesh_provider_->load_surface_from_file("/Users/kraemer/Data/surface/grid_tri_high.off");

		simulation::shallow_water_get_attributes(*domain2D_, sw_attributes);
		simulation::shallow_water_init_attributes(*domain2D_, sw_attributes, sw_context);

		vertex_water_position_ = add_attribute<Vec3, Vertex>(*domain2D_, "water_position");
		vertex_water_position_->copy(sw_attributes.vertex_position_.get());

		update_render_data();

		running_ = false;

		timer_connection_ = boost::synapse::connect<App::timer_tick>(&app_, [this] ()
		{
			update_render_data();
			for (View* v : linked_views_)
				v->request_update();
		});
	}

	void start()
	{
		running_ = true;

		launch_thread([this] ()
		{
			while (this->running_)
			{
				shallow_water_execute_time_step(*domain2D_, sw_attributes, sw_context);
				if (sw_context.t_ == sw_context.t_max_)
					stop();
			}
		});

		app_.start_timer(50, [this] () -> bool { return !running_; });
	}

	void stop()
	{
		running_ = false;
	}

	void update_render_data()
	{
		foreach_cell(*domain2D_, [&] (Vertex v) -> bool
		{
			Scalar h = 0.0;
			uint32 nbf = 0;
			foreach_incident_face(*domain2D_, v, [&] (Face f) -> bool
			{
				h += value<Scalar>(*domain2D_, sw_attributes.face_h_, f);
				++nbf;
				return true;
			});
			value<Vec3>(*domain2D_, vertex_water_position_, v)[2] = h / nbf;
			return true;
		});
		
		mesh_provider_->emit_attribute_changed(domain2D_, vertex_water_position_.get());
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
			{
				shallow_water_execute_time_step(*domain2D_, sw_attributes, sw_context);
				update_render_data();
				for (View* v : linked_views_)
					v->request_update();
			}
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

	std::shared_ptr<Attribute<Vec3>> vertex_water_position_;

	simulation::ShallowWaterAttributes<MESH> sw_attributes;
	simulation::ShallowWaterContext sw_context;

	std::shared_ptr<boost::synapse::connection> timer_connection_;
	bool running_ = false;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SHALLOW_WATER_H_

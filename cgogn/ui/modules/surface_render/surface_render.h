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

#ifndef CGOGN_MODULE_SURFACE_RENDER_H_
#define CGOGN_MODULE_SURFACE_RENDER_H_

#include <cgogn/ui/modules/surface_render/cgogn_module_surface_render_export.h>

#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/mesh_render.h>
#include <cgogn/rendering/vbo_update.h>

#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_phong.h>
#include <cgogn/rendering/shaders/shader_vector_per_vertex.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>

namespace cgogn
{

namespace ui
{

class App;
class View;

class CGOGN_MODULE_SURFACE_RENDER_EXPORT SurfaceRender : public Module
{
    using Mesh = CMap2;

    template <typename T>
    using Attribute = typename mesh_traits<Mesh>::Attribute<T>;

    using Vertex = typename mesh_traits<Mesh>::Vertex;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

public:

	SurfaceRender(const App& app);
	~SurfaceRender();

	void init();

	void update_data(const Mesh& m);

protected:

	void draw(View* view) override;
    void interface() override;

private:

	struct Parameters
	{
		Parameters() : mesh_(nullptr)
		{}

		Parameters(const Mesh* m, MeshData* m_data) :
			mesh_(m),
			mesh_data_(m_data),
			initialized_(false),
			render_vertices_(false),
			render_edges_(false),
			render_faces_(true),
			phong_shading_(false),
			vertex_scale_factor_(1.0)
		{
			param_point_sprite_ = rendering::ShaderPointSprite::generate_param();
			// param_point_sprite_->set_vbos(vbo_position_.get());
			param_point_sprite_->color_ = rendering::GLColor(1, 0, 0, 1);

			param_edge_ = rendering::ShaderBoldLine::generate_param();
			// param_edge_->set_vbos(vbo_position_.get());
			param_edge_->color_ = rendering::GLColor(1, 1, 0, 1);
			param_edge_->width_= 2.5f;

			param_flat_ =  rendering::ShaderFlat::generate_param();
			// param_flat_->set_vbos(vbo_position_.get());
			param_flat_->front_color_ = rendering::GLColor(0, 0.8f, 0, 1);
			param_flat_->back_color_ = rendering::GLColor(0, 0, 0.8f, 1);
			param_flat_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);

			param_phong_ = rendering::ShaderPhong::generate_param();
			// param_phong_->set_vbos(vbo_position_.get(), vbo_normal_.get());
			param_phong_->front_color_ = rendering::GLColor(0, 0.8f, 0, 1);
			param_phong_->back_color_ = rendering::GLColor(0, 0, 0.8f, 1);
			param_phong_->ambiant_color_ = rendering::GLColor(0.1f, 0.1f, 0.1f, 1);
			param_phong_->specular_coef_ = 250.0f;
		}

		const Mesh* mesh_;
		MeshData* mesh_data_;
		bool initialized_;

		std::string vertex_position_name_;
		std::string vertex_normal_name_;

		std::unique_ptr<rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<rendering::ShaderFlat::Param> param_flat_;
		std::unique_ptr<rendering::ShaderPhong::Param> param_phong_;

		bool render_vertices_;
		bool render_edges_;
		bool render_faces_;
		bool phong_shading_;
		
		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
	};

	const Mesh* selected_mesh_;
	std::unordered_map<const Mesh*, Parameters> parameters_;
	ui::MeshProvider* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_RENDER_H_

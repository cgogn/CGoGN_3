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
class CMapProvider;

class CGOGN_MODULE_SURFACE_RENDER_EXPORT SurfaceRender : public Module
{
    using Mesh = CMap2;

    template <typename T>
    using AttributePtr = typename cgogn::mesh_traits<Mesh>::AttributePtr<T>;
    using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;

    using Vec3 = cgogn::geometry::Vec3;
    using Scalar = cgogn::geometry::Scalar;

public:

	SurfaceRender(const App& app);
	~SurfaceRender();

	void init();

	void update(Mesh* m, const AttributePtr<Vec3>& vertex_position,  const AttributePtr<Vec3>& vertex_normal);

protected:

	void draw(View* view) override;
    void interface() override;

private:

	struct Parameters
	{
		Parameters() : mesh_(nullptr)
		{}

		Parameters(Mesh* m) :
			mesh_(m),
			initialized_(false),
			render_vertices_(false),
			render_edges_(false),
			render_faces_(true),
			phong_shading_(false),
			vertex_scale_factor_(1.0)
		{
			render_ = std::make_unique<cgogn::rendering::MeshRender>();

			vbo_position_ = std::make_unique<cgogn::rendering::VBO>();
			vbo_normal_ = std::make_unique<cgogn::rendering::VBO>();

			param_point_sprite_ = cgogn::rendering::ShaderPointSprite::generate_param();
			param_point_sprite_->set_vbos(vbo_position_.get());
			param_point_sprite_->color_ = cgogn::rendering::GLColor(1, 0, 0, 1);

			param_edge_ = cgogn::rendering::ShaderBoldLine::generate_param();
			param_edge_->set_vbos(vbo_position_.get());
			param_edge_->color_ = cgogn::rendering::GLColor(1, 1, 0, 1);
			param_edge_->width_= 2.5f;

			param_flat_ =  cgogn::rendering::ShaderFlat::generate_param();
			param_flat_->set_vbos(vbo_position_.get());
			param_flat_->front_color_ = cgogn::rendering::GLColor(0, 0.8f, 0, 1);
			param_flat_->back_color_ = cgogn::rendering::GLColor(0, 0, 0.8f, 1);
			param_flat_->ambiant_color_ = cgogn::rendering::GLColor(0.1f, 0.1f, 0.1f, 1);

			param_phong_ = cgogn::rendering::ShaderPhong::generate_param();
			param_phong_->set_vbos(vbo_position_.get(), vbo_normal_.get());
			param_phong_->front_color_ = cgogn::rendering::GLColor(0, 0.8f, 0, 1);
			param_phong_->back_color_ = cgogn::rendering::GLColor(0, 0, 0.8f, 1);
			param_phong_->ambiant_color_ = cgogn::rendering::GLColor(0.1f, 0.1f, 0.1f, 1);
			param_phong_->specular_coef_ = 250.0f;
		}

		Mesh* mesh_;
		bool initialized_;

		std::unique_ptr<cgogn::rendering::MeshRender> render_;
		Vec3 bb_min_, bb_max_;

		std::unique_ptr<cgogn::rendering::VBO> vbo_position_;
		std::unique_ptr<cgogn::rendering::VBO> vbo_normal_;

		std::unique_ptr<cgogn::rendering::ShaderPointSprite::Param> param_point_sprite_;
		std::unique_ptr<cgogn::rendering::ShaderBoldLine::Param> param_edge_;
		std::unique_ptr<cgogn::rendering::ShaderFlat::Param> param_flat_;
		std::unique_ptr<cgogn::rendering::ShaderPhong::Param> param_phong_;

		bool render_vertices_;
		bool render_edges_;
		bool render_faces_;
		bool phong_shading_;
		
		float32 vertex_scale_factor_;
		float32 vertex_base_size_;
	};

	Mesh* selected_mesh_;
	AttributePtr<Vec3> selected_vertex_position_;
	AttributePtr<Vec3> selected_vertex_normal_;
	std::unordered_map<Mesh*, Parameters> parameters_;
	cgogn::ui::CMapProvider* cmap_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_RENDER_H_

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

#ifndef CGOGN_RENDERING_SHADER_TRANSP_VOLUMES_H_
#define CGOGN_RENDERING_SHADER_TRANSP_VOLUMES_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/shaders/vbo.h>

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/ear_triangulation.h>
#include <cgogn/geometry/types/geometry_traits.h>

#include <GLColor>
#include <QOpenGLFramebufferObject>
#include <QOpenGLFunctions_3_3_Core>

namespace cgogn
{

namespace rendering
{

class ShaderTransparentVolumes;

class ShaderParamTransparentVolumes : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	using ShaderType = ShaderTransparentVolumes;

	GLColor color_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;
	QVector3D light_position_;
	float32 explode_factor_;

	bool bf_culling_;
	bool lighted_;

	ShaderParamTransparentVolumes(ShaderTransparentVolumes* sh);

	void set_position_vbo(VBO* vbo_pos);
};

class CGOGN_RENDERING_EXPORT ShaderTransparentVolumes : public ShaderProgram
{
	friend class ShaderParamTransparentVolumes;

protected:
	static const char* vertex_shader_source_;
	static const char* geometry_shader_source_;
	static const char* fragment_shader_source_;

	// uniform ids
	GLint unif_expl_v_;
	GLint unif_light_position_;
	GLint unif_plane_clip_;
	GLint unif_plane_clip2_;
	GLint unif_color_;
	GLint unif_bf_culling_;
	GLint unif_lighted_;
	GLint unif_layer_;
	GLint unif_depth_texture_sampler_;
	GLint unif_rgba_texture_sampler_;

public:
	using Self = ShaderTransparentVolumes;
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ShaderTransparentVolumes);

	enum
	{
		ATTRIB_POS = 0,
	};

	void set_explode_volume(float32 x);
	void set_light_position(const QVector3D& l);
	void set_plane_clip(const GLVec4& plane);
	void set_plane_clip2(const GLVec4& plane);
	void set_color(const GLColor& rgb);

	void set_bf_culling(bool cull);
	void set_lighted(bool lighted);
	void set_layer(int layer);
	void set_rgba_sampler(GLuint rgba_samp);
	void set_depth_sampler(GLuint depth_samp);

	using Param = ShaderParamTransparentVolumes;

	static std::unique_ptr<Param> generate_param();

	static ShaderTransparentVolumes* get_instance();

protected:
	ShaderTransparentVolumes();

	static ShaderTransparentVolumes* instance_;
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADER_TRANSP_VOLUMES_H_

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

#ifndef CGOGN_RENDERING_SHADER_TRANSP_PHONG_H_
#define CGOGN_RENDERING_SHADER_TRANSP_PHONG_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/shaders/vbo.h>

#include <GLColor>
#include <QOpenGLFramebufferObject>
#include <QOpenGLFunctions>

namespace cgogn
{

namespace rendering
{

// forward
class ShaderPhongTransp;

class CGOGN_RENDERING_EXPORT ShaderParamPhongTransp : public ShaderParam
{
protected:
	void set_uniforms() override;

public:
	using ShaderType = ShaderPhongTransp;

	GLColor front_color_;
	GLColor back_color_;
	GLColor ambiant_color_;
	GLColor specular_color_;
	GLfloat specular_coef_;
	QVector3D light_pos_;
	bool bf_culling_;

	ShaderParamPhongTransp(ShaderPhongTransp* sh);

	void set_position_vbo(VBO* vbo_pos);

	void set_normal_vbo(VBO* vbo_normal);

	void set_alpha(int alpha);
};

class CGOGN_RENDERING_EXPORT ShaderPhongTransp : public ShaderProgram
{
	friend class ShaderParamPhongTransp;

protected:
	static const char* vertex_shader_source_;
	static const char* fragment_shader_source_;

	// uniform ids
	GLint unif_front_color_;
	GLint unif_back_color_;
	GLint unif_ambiant_color_;
	GLint unif_specular_color_;
	GLint unif_specular_coef_;
	GLint unif_light_position_;
	GLint unif_bf_culling_;
	GLint unif_layer_;
	GLint unif_depth_texture_sampler_;
	GLint unif_rgba_texture_sampler_;

public:
	using Self = ShaderPhongTransp;
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ShaderPhongTransp);

	enum
	{
		ATTRIB_POS = 0,
		ATTRIB_NORM
	};

	/**
	 * @brief set current front color
	 * @param rgb
	 */
	void set_front_color(const GLColor& rgb);

	/**
	 * @brief set current front color
	 * @param rgb
	 */
	void set_back_color(const GLColor& rgb);

	/**
	 * @brief set current ambiant color
	 * @param rgb
	 */
	void set_ambiant_color(const GLColor& rgb);

	/**
	 * @brief set current specular color
	 * @param rgb
	 */
	void set_specular_color(const GLColor& rgb);

	/**
	 * @brief set current specular coefficient
	 * @param rgb
	 */
	void set_specular_coef(float32 coef);

	/**
	 * @brief set light position relative to screen
	 * @param l light position
	 */
	void set_light_position(const QVector3D& l);

	void set_bf_culling(bool cull);

	void set_layer(int layer);

	void set_rgba_sampler(GLuint rgba_samp);

	void set_depth_sampler(GLuint depth_samp);

	using Param = ShaderParamPhongTransp;

	static std::unique_ptr<Param> generate_param();

	static ShaderPhongTransp* get_instance();

private:
	ShaderPhongTransp();
	static ShaderPhongTransp* instance_;
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADER_TRANSP_PHONG_H_

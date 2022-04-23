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

#ifndef CGOGN_RENDERING_SHAPE3D_DRAWER_H_
#define CGOGN_RENDERING_SHAPE3D_DRAWER_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shader_program.h>

namespace cgogn
{

namespace rendering
{

DECLARE_SHADER_CLASS(Cylinder, false, CGOGN_STR(Cylinder))

class CGOGN_RENDERING_EXPORT ShaderParamCylinder : public ShaderParam
{
	void set_uniforms() override;

public:
	int32 nb_;
	GLColor color_;
	float32 roughness_;
	float32 shininess_;
	GLVec3 light_position_;

	using ShaderType = ShaderCylinder;

	ShaderParamCylinder(ShaderType* sh);

	~ShaderParamCylinder() override;
};

DECLARE_SHADER_CLASS(Sphere, false, CGOGN_STR(Sphere))

class CGOGN_RENDERING_EXPORT ShaderParamSphere : public ShaderParam
{
	void set_uniforms() override;

public:
	int32 nb_;
	GLColor color_;
	float32 roughness_;
	float32 shininess_;
	GLVec3 light_position_;

	using ShaderType = ShaderSphere;

	ShaderParamSphere(ShaderType* sh);

	 ~ShaderParamSphere() override;
};


DECLARE_SHADER_CLASS(Cone, false, CGOGN_STR(Cone))

class CGOGN_RENDERING_EXPORT ShaderParamCone : public ShaderParam
{
	void set_uniforms() override;

public:
	int32 nb_;
	GLColor color_;
	float32 roughness_;
	float32 shininess_;
	GLVec3 light_position_;

	using ShaderType = ShaderCone;

	ShaderParamCone(ShaderType* sh);

	 ~ShaderParamCone() override;
};


class Shape3DDrawer
{
	int32 nb_subd_;
	std::unique_ptr<ShaderCylinder::Param> param_cylinder_;
	std::unique_ptr<ShaderSphere::Param> param_sphere_;
	std::unique_ptr<ShaderCone::Param> param_cone_;

	static Shape3DDrawer* instance_;

public:
	Shape3DDrawer();
	~Shape3DDrawer();

	static Shape3DDrawer* instance();

	void update_material_cylinder(const GLColor& col, float32 roughness, float32 shininess);

	void update_material_sphere(const GLColor& col, float32 roughness, float32 shininess);

	void update_material_cone(const GLColor& col, float32 roughness, float32 shininess);

	void update_light_position(const GLVec3& lp);

	void draw_cylinder(const GLMat4& projection, const GLMat4& view, const GLMat4& transfo);
	void draw_sphere(const GLMat4& projection, const GLMat4& view, const GLMat4& transfo);
	void draw_cone(const GLMat4& projection, const GLMat4& view, const GLMat4& transfo);
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHAPE3_DRAWER_H_




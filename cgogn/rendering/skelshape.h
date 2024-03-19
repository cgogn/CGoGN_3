<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<</*******************************************************************************
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

#ifndef CGOGN_RENDERING_SKEL_SHAPE_DRAWER_H_
#define CGOGN_RENDERING_SKEL_SHAPE_DRAWER_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shader_program.h>

namespace cgogn
{

namespace rendering
{
/*
ajouter:
void ShaderParam::set_vbos_inst(const std::vector<VBO*>& vbos)
void VBO::associate(GLuint attrib, int32 stride = 0, uint32 first = 0, uint32 divisor = 0 )
*/	
class ShaderParamSkelShape : public ShaderParam
{
public:
	inline virtual void set_subdiv(int32 nb){};
	inline virtual void draw(const GLMat4& projection, const GLMat4& view) = 0;

	int nb_;
	GLColor color_;
	float roughness_;
	float shininess_;
	GLVec3 light_position_;

	std::unique_ptr<rendering::VAO> inst_vao_;

	ShaderParamSkelShape(ShaderProgram* prg);

	virtual void draw_inst(const GLMat4& projection, const GLMat4& view) = 0;

};


DECLARE_SHADER_CLASS(SkelCone, false, CGOGN_STR(SkelCone))

class CGOGN_RENDERING_EXPORT ShaderParamSkelSphere : public ShaderParamSkelShape
{

	void set_uniforms() override;

public:

	using ShaderType = ShaderSphere;

	inline ShaderParamSkelSphere(ShaderType* sh):
	ShaderParamSkelShape(sh) 
	{}

	inline ~ShaderParamSkelSphere() override {}

	void draw_inst(const GLMat4& projection, const GLMat4& view) override;
};

class SkelSphereDrawer
{
    std::unique_ptr<rendering::VBO> vbo_;
	std::unique_ptr<ShaderSkelSphere::Param> param_sphere_;
public:
	inline void set_spheres(const std::vector<GLVec4>& P1)
    {
       update_vbo(P1,vbo_.get());
       param_sphere_->set_vbo({vbo_.get()});
    }

	void draw(const GLMat4& projection, const GLMat4& view) override;
};

// DECLARE_SHADER_CLASS(SkelCone, false, CGOGN_STR(SkelCone))

// class CGOGN_RENDERING_EXPORT ShaderParamSkelCone : public ShaderParamSkelShape
// {

// 	void set_uniforms() override;

// public:

// 	using ShaderType = ShaderCone;

// 	inline ShaderParamSkelCone(ShaderType* sh):
// 	ShaderParamSkelShape(sh) 
// 	{}

// 	inline ~ShaderParamSkelCone() override {}

// 	void draw(const GLMat4& projection, const GLMat4& view) override;
// };




// class SkelConeDrawer
// {
// 	std::unique_ptr<ShaderSkelCone::Param> param_cone_;
//     rendering::VBO> vbo1_;
//     rendering::VBO> vbo2_;
// public:
//     inline void set_cones(const std::vector<GLVec4>& P1, const std::vector<GLVec4>& P2)
//     {
//        update_vbo(P1,&vbo_)
//        set_vbo({&vbo_});
//     }
// };

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHAPE3_DRAWER_H_





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

	
class ShaderParamShape : public ShaderParam
{
public:
	inline virtual void set_subdiv(int32 nb){};
	inline virtual void draw(const GLMat4& projection, const GLMat4& view) = 0;

	GLColor color_;
	float roughness_;
	float shininess_;
	GLVec3 light_position_;

	ShaderParamShape(ShaderProgram* prg);

};

class ShaderParamShapeSub : public ShaderParamShape
{
public:
	int nb_;

	inline void set_subdiv(int32 nb)
	{
		nb_ = nb;
	}

	ShaderParamShapeSub(ShaderProgram* prg);

};


DECLARE_SHADER_CLASS(Cylinder, false, CGOGN_STR(Cylinder))

class CGOGN_RENDERING_EXPORT ShaderParamCylinder : public ShaderParamShapeSub
{
	void set_uniforms() override;

public:

	using ShaderType = ShaderCylinder;

	ShaderParamCylinder(ShaderType* sh);

	~ShaderParamCylinder() override;

	void draw(const GLMat4& projection, const GLMat4& view) override;
};


DECLARE_SHADER_CLASS(Sphere, false, CGOGN_STR(Sphere))

class CGOGN_RENDERING_EXPORT ShaderParamSphere : public ShaderParamShapeSub
{
	void set_uniforms() override;

public:

	using ShaderType = ShaderSphere;

	ShaderParamSphere(ShaderType* sh);

	~ShaderParamSphere() override;

	void draw(const GLMat4& projection, const GLMat4& view) override;
};


DECLARE_SHADER_CLASS(Cone, false, CGOGN_STR(Cone))

class CGOGN_RENDERING_EXPORT ShaderParamCone : public ShaderParamShapeSub
{
	void set_uniforms() override;

public:

	using ShaderType = ShaderCone;

	ShaderParamCone(ShaderType* sh);

	~ShaderParamCone() override;

	void draw(const GLMat4& projection, const GLMat4& view) override;
};


DECLARE_SHADER_CLASS(Cube, false, CGOGN_STR(Cube))

class CGOGN_RENDERING_EXPORT ShaderParamCube : public ShaderParamShape
{
	void set_uniforms() override;
public:

	using ShaderType = ShaderCube;

	ShaderParamCube(ShaderType* sh);

	~ShaderParamCube() override;

	void draw(const GLMat4& projection, const GLMat4& view) override;
};


class ShapeDrawer
{
public:
	enum SHAPE
	{
		CYLINDER=0,
		SPHERE,
		CONE,
		CUBE,
		NBSHAPES
	};

	static ShapeDrawer* instance();

protected:
	int32 nb_subd_;

	std::array<std::unique_ptr<ShaderParamShape>, SHAPE::NBSHAPES> param_;

	static ShapeDrawer* instance_;

	ShapeDrawer();

	GLMat4 points_2_transfo(const GLVec3& p1, const GLVec3& p2, float32 radius);

public:


	~ShapeDrawer();

	/**
	* @brief get a reference to color of shape s
	* @param s shape [CYLINDER,SPHERE,CONE,CUBE]
	* @return the refence
	*/
	inline GLColor& color(SHAPE s)
	{
		return param_[s]->color_;
	}

	/**
	* @brief get a reference to roughness phong exponent [1-250] of shape s
	* @param s shape [CYLINDER,SPHERE,CONE,CUBE]
	* @return the refence
	*/
	inline float& roughness(SHAPE s)
	{
		return param_[s]->roughness_;
	}

	/**
	* @brief get a reference to shiness [0,1] of shape s
	* @param s shape [CYLINDER,SPHERE,CONE,CUBE]
	* @return the refence
	*/
	inline float& shininess(SHAPE s)
	{
		return param_[s]->shininess_;
	}

	/**
	* @brief update the material of all shapes
	* @param colo RGBA color
	* @param roughness roughness phong exponent [1-250]
	* @param shininess shininess factor [0,1]
	*/
	void update_material(const GLColor& col, float32 roughness, float32 shininess);

	/**
	* @brief update the light position for all shapes
	* @param lp
	*/
	void update_light_position(const GLVec3& lp);

	/**
	 * @brief update the subdivion for all smooth shapes
	 * @param nb
	 */
	void update_subdivision(int32 nb);


	/**
	* @brief draw
	* @param s shape [CYLINDER,SPHERE,CONE,CUBE]
	* @param projection matrix from View
	* @param model_view matrix from View multiplied by local transfo
	*/
	void draw(SHAPE s, const GLMat4& projection, const GLMat4& model_view);

	/**
	* @brief draw a sphere at position p of radius size
	* @param projection matrix from View
	* @param view matrix from View
	* @param p position
	* @param radius
	*/
	void drawSphere(const GLMat4& projection, const GLMat4& view, const GLVec3& p, float32 radius);

	/**
	* @brief draw a cylinder with given radius from point p1 to point p2
	* @param projection matrix from View
	* @param view matrix from View
	* @param p1 start point
	* @param p2 end point
	* @param radius
	*/
	void drawCylinder(const GLMat4& projection, const GLMat4& view, const GLVec3& p1, const GLVec3& p2, float32 radius);

	/**
	* @brief draw a cone with given radius from point p1 to point p2
	* @param projection matrix from View
	* @param view matrix from View
	* @param p1 start point
	* @param p2 end point
	* @param radius
	*/
	void drawCone(const GLMat4& projection, const GLMat4& view, const GLVec3& p1, const GLVec3& p2, float32 radius);
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHAPE3_DRAWER_H_




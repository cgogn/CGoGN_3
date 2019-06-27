/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#ifndef CGOGN_RENDERING_SHADERS_SCALARPERVERTEX_H_
#define CGOGN_RENDERING_SHADERS_SCALARPERVERTEX_H_

#include <cgogn/rendering_pureGL/cgogn_rendering_puregl_export.h>
#include <cgogn/rendering_pureGL/shaders/shader_program.h>

namespace cgogn
{

namespace rendering_pgl
{

// forward
class ShaderParamScalarPerVertex;

class CGOGN_RENDERING_PUREGL_EXPORT ShaderScalarPerVertex : public ShaderProgram
{
public:
	using  Self  = ShaderScalarPerVertex;
	using  Param = ShaderParamScalarPerVertex;
	friend Param;

	enum ColorMap
	{
		BWR = 0,
		CWR,
		BCGYR,
		BGR
	};
protected:
	ShaderScalarPerVertex();
	CGOGN_NOT_COPYABLE_NOR_MOVABLE(ShaderScalarPerVertex);
	static Self* instance_;
public:
	inline static std::unique_ptr<Param> generate_param()
	{
		if (!instance_)
		{
			instance_ = new Self();
			ShaderProgram::register_instance(instance_);
		}
		return cgogn::make_unique<Param>(instance_);
	}
};


class CGOGN_RENDERING_PUREGL_EXPORT ShaderParamScalarPerVertex : public ShaderParam
{
	inline void set_uniforms() override
	{
		shader_->set_uniforms_values(
					color_map_,
					expansion_,
					min_value_,
					max_value_,
					show_iso_lines_,
					nb_iso_levels_);
	}

public:
	ShaderScalarPerVertex::ColorMap color_map_;
	int32 expansion_;
	float32 min_value_;
	float32 max_value_;
	bool show_iso_lines_;
	int32 nb_iso_levels_;

	using LocalShader = ShaderScalarPerVertex;

	ShaderParamScalarPerVertex(LocalShader* sh) :
		ShaderParam(sh),
		color_map_(ShaderScalarPerVertex::BWR),
		expansion_(0),
		min_value_(.0f),
		max_value_(1.0f),
		show_iso_lines_(false),
		nb_iso_levels_(10)
	{}

	inline ~ShaderParamScalarPerVertex() override {}

	inline void set_vbos(VBO* vbo_pos, VBO* vbo_scalar)
	{
		bind_vao();
		associate_vbos(vbo_pos,vbo_scalar);
		release_vao();
	}

};


} // namespace rendering_pgl
} // namespace cgogn
#endif // CGOGN_RENDERING_SHADERS_SCALARPERVERTEX_H_

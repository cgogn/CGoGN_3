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

#include <cgogn/core/utils/unique_ptr.h>
#include <cgogn/rendering/shader_program.h>

namespace cgogn
{

namespace rendering
{

ShaderParam::ShaderParam(ShaderProgram* prg) :
	shader_(prg)
{
	vao_ = cgogn::make_unique<QOpenGLVertexArrayObject>();
	vao_->create();
}

ShaderParam::~ShaderParam()
{}

void ShaderParam::bind_vao_only(bool with_uniforms)
{
	if (with_uniforms)
		set_uniforms();
	vao_->bind();
}

void ShaderParam::release_vao_only()
{
	vao_->release();
}

void ShaderParam::bind(const QMatrix4x4& proj, const QMatrix4x4& mv)
{
	shader_->bind();
	shader_->set_matrices(proj,mv);
	set_uniforms();
	vao_->bind();
}

void ShaderParam::release()
{
	vao_->release();
	shader_->release();
}

/*****************************************************************************/

std::vector<ShaderProgram*>* ShaderProgram::instances_ = nullptr;

void ShaderProgram::register_instance(ShaderProgram* sh)
{
	if (instances_ == nullptr)
	{
		instances_ = new std::vector<ShaderProgram*>;
		instances_->reserve(256);
	}
	
	auto it = std::find(instances_->begin(), instances_->end(), sh);
	if (it == instances_->end())
		instances_->push_back(sh);
}

void ShaderProgram::clean_all()
{
	if (instances_ != nullptr)
	{
		for (auto* ptr : *instances_)
			delete ptr;
		delete instances_;
		instances_ = nullptr;
	}
}

ShaderProgram::~ShaderProgram()
{}

void ShaderProgram::get_matrices_uniforms()
{
	unif_mv_matrix_ = prg_.uniformLocation("model_view_matrix");
	unif_projection_matrix_ = prg_.uniformLocation("projection_matrix");
	unif_normal_matrix_ = prg_.uniformLocation("normal_matrix");
}

void ShaderProgram::set_matrices(const QMatrix4x4& proj, const QMatrix4x4& mv)
{
	if (unif_projection_matrix_ >= 0)
		prg_.setUniformValue(unif_projection_matrix_, proj);
	if (unif_mv_matrix_ >= 0)
		prg_.setUniformValue(unif_mv_matrix_, mv);

	if (unif_normal_matrix_ >= 0)
	{
		QMatrix3x3 normal_matrix = mv.normalMatrix();
		prg_.setUniformValue(unif_normal_matrix_, normal_matrix);
	}
}

void ShaderProgram::set_view_matrix(const QMatrix4x4& mv)
{
	prg_.setUniformValue(unif_mv_matrix_, mv);

	if (unif_normal_matrix_ >= 0)
	{
		QMatrix3x3 normal_matrix = mv.normalMatrix();
		prg_.setUniformValue(unif_normal_matrix_, normal_matrix);
	}
}

} // namespace rendering

} // namespace cgogn

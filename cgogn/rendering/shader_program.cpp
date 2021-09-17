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

#include <cgogn/rendering/shader_program.h>

namespace cgogn
{

namespace rendering
{

/*****************************************************************************/
// Shader
/*****************************************************************************/

void Shader::compile(const std::string& src, const std::string& prg_name)
{
	const char* csrc = src.c_str();
	glShaderSource(id_, 1, &csrc, nullptr);
	glCompileShader(id_);

	int infologLength = 0;
	int charsWritten = 0;
	char* infoLog;

	glGetShaderiv(id_, GL_INFO_LOG_LENGTH, &infologLength);

	if (infologLength > 1)
	{
		infoLog = (char*)malloc(infologLength + 1);
		glGetShaderInfoLog(id_, infologLength, &charsWritten, infoLog);

		std::cerr << "----------------------------------------" << std::endl
				  << "compilation de " << prg_name << " : " << std::endl
				  << infoLog << std::endl
				  << "--------" << std::endl;

		std::string errors(infoLog);
		std::istringstream sserr(errors);
		std::vector<int> error_lines;
		std::string line;
		std::getline(sserr, line);
		while (!sserr.eof())
		{
			std::size_t a = 0;
			while ((a < uint32(line.size())) && (line[a] >= '0') && (line[a] <= '9'))
				a++;
			std::size_t b = a + 1;
			while ((b < uint32(line.size())) && (line[b] >= '0') && (line[b] <= '9'))
				b++;
			if (b < uint32(line.size()))
			{
				int ln = std::stoi(line.substr(a + 1, b - a - 1));
				error_lines.push_back(ln);
			}
			std::getline(sserr, line);
		}

		free(infoLog);

		char* source = new char[16 * 1024];
		GLsizei length;
		glGetShaderSource(id_, 16 * 1024, &length, source);
		std::string sh_src(source);
		std::istringstream sssrc(sh_src);
		int l = 1;
		while (!sssrc.eof())
		{
			std::getline(sssrc, line);
			std::cerr.width(3);
			auto it = std::find(error_lines.begin(), error_lines.end(), l);
			if (it != error_lines.end())
				std::cerr << "\033[41m\033[37m"
						  << "EEEEEE" << line << "\033[m" << std::endl;
			else
				std::cerr << l << " : " << line << std::endl;
			l++;
		}
		std::cerr << "----------------------------------------" << std::endl;
	}
}

/*****************************************************************************/
// ShaderProgram
/*****************************************************************************/

std::vector<ShaderProgram*>* ShaderProgram::instances_ = nullptr;

ShaderProgram::ShaderProgram()
	: vertex_shader_(nullptr), fragment_shader_(nullptr), geometry_shader_(nullptr), nb_attributes_(0)
{
	id_ = glCreateProgram();
}

ShaderProgram::~ShaderProgram()
{
	if (vertex_shader_)
		delete vertex_shader_;
	if (geometry_shader_)
		delete geometry_shader_;
	if (fragment_shader_)
		delete fragment_shader_;

	glDeleteProgram(id_);
}

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

void ShaderProgram::get_matrices_uniforms()
{
	uniform_mvp_matrix_ = glGetUniformLocation(id_, "mvp_matrix");
	uniform_mv_matrix_ = glGetUniformLocation(id_, "model_view_matrix");
	uniform_projection_matrix_ = glGetUniformLocation(id_, "projection_matrix");
	uniform_normal_matrix_ = glGetUniformLocation(id_, "normal_matrix");
}

void ShaderProgram::set_matrices(const GLMat4d& projection, const GLMat4d& mv)
{
	if (uniform_mvp_matrix_ >= 0)
	{
		GLMat4d mvp = projection * mv;
		GLMat4 m = mvp.cast<float>();
		glUniformMatrix4fv(uniform_mvp_matrix_, 1, false, m.data());
	}
	if (uniform_projection_matrix_ >= 0)
	{
		GLMat4 m = projection.cast<float>();
		glUniformMatrix4fv(uniform_projection_matrix_, 1, false, m.data());
	}
	if (uniform_mv_matrix_ >= 0)
	{
		GLMat4 m = mv.cast<float>();
		glUniformMatrix4fv(uniform_mv_matrix_, 1, false, m.data());
	}
	if (uniform_normal_matrix_ >= 0)
	{
		Eigen::Affine3d t(mv);
		GLMat3 normal_matrix = t.linear().inverse().transpose().matrix().cast<float>();
		glUniformMatrix3fv(uniform_normal_matrix_, 1, false, normal_matrix.data());
	}
}

void ShaderProgram::set_matrices(const GLMat4& projection, const GLMat4& mv)
{
	if (uniform_mvp_matrix_ >= 0)
	{
		GLMat4 m = projection * mv;
		glUniformMatrix4fv(uniform_mvp_matrix_, 1, false, m.data());
	}
	if (uniform_projection_matrix_ >= 0)
		glUniformMatrix4fv(uniform_projection_matrix_, 1, false, projection.data());
	if (uniform_mv_matrix_ >= 0)
		glUniformMatrix4fv(uniform_mv_matrix_, 1, false, mv.data());
	if (uniform_normal_matrix_ >= 0)
	{
		Eigen::Affine3d t(mv.cast<float64>());
		GLMat3 normal_matrix = t.linear().matrix().inverse().transpose().cast<float32>();
		glUniformMatrix3fv(uniform_normal_matrix_, 1, false, normal_matrix.data());
	}
}

void ShaderProgram::load(const std::string& vert_src, const std::string& frag_src)
{
	vertex_shader_ = new Shader(GL_VERTEX_SHADER);
	vertex_shader_->compile(vert_src, name());

	fragment_shader_ = new Shader(GL_FRAGMENT_SHADER);
	fragment_shader_->compile(frag_src, name());

	glAttachShader(id_, vertex_shader_->id());
	glAttachShader(id_, fragment_shader_->id());

	glLinkProgram(id_);

	// puis detache (?)
	glDetachShader(id_, fragment_shader_->id());
	glDetachShader(id_, vertex_shader_->id());

	// glValidateProgram(id_);
	// // Print log if needed
	// GLint infologLength = 0;
	// glGetProgramiv(id_, GL_LINK_STATUS, &infologLength);
	// if (infologLength != GL_TRUE)
	// 	std::cerr << "PB GL_LINK_STATUS load " << name() << " " << infologLength << std::endl;
	// glGetProgramiv(id_, GL_VALIDATE_STATUS, &infologLength);
	// if (infologLength != GL_TRUE)
	// 	std::cerr << "PB GL_VALIDATE_STATUS load " << name() << " " << infologLength << std::endl;
	// glGetProgramiv(id_, GL_INFO_LOG_LENGTH, &infologLength);
	// if (infologLength > 1)
	// {
	// 	char* infoLog = new char[infologLength];
	// 	int charsWritten = 0;
	// 	glGetProgramInfoLog(id_, infologLength, &charsWritten, infoLog);
	// 	std::cerr << "Link message: " << infoLog << std::endl;
	// 	delete[] infoLog;
	// }

	get_matrices_uniforms();
}

void ShaderProgram::load3(const std::string& vert_src, const std::string& frag_src, const std::string& geom_src)
{
	vertex_shader_ = new Shader(GL_VERTEX_SHADER);
	vertex_shader_->compile(vert_src, name());

	geometry_shader_ = new Shader(GL_GEOMETRY_SHADER);
	geometry_shader_->compile(geom_src, name());

	fragment_shader_ = new Shader(GL_FRAGMENT_SHADER);
	fragment_shader_->compile(frag_src, name());

	glAttachShader(id_, vertex_shader_->id());
	glAttachShader(id_, geometry_shader_->id());
	glAttachShader(id_, fragment_shader_->id());

	glLinkProgram(id_);

	// puis detache (?)
	glDetachShader(id_, fragment_shader_->id());
	glDetachShader(id_, geometry_shader_->id());
	glDetachShader(id_, vertex_shader_->id());

	// glValidateProgram(id_);
	// // Print log if needed
	// GLint infologLength = 0;
	// glGetProgramiv(id_, GL_LINK_STATUS, &infologLength);
	// if (infologLength != GL_TRUE)
	// 	std::cerr << "PB GL_LINK_STATUS load " << name() << " " << infologLength << std::endl;
	// glGetProgramiv(id_, GL_VALIDATE_STATUS, &infologLength);
	// if (infologLength != GL_TRUE)
	// 	std::cerr << "PB GL_VALIDATE_STATUS load " << name() << " " << infologLength << std::endl;
	// glGetProgramiv(id_, GL_ATTACHED_SHADERS, &infologLength);
	// std::cerr << "NB ATTACHED SHADERS in " << name() << " " << infologLength << std::endl;
	// glGetProgramiv(id_, GL_ACTIVE_UNIFORMS, &infologLength);
	// std::cerr << "NB ACTIVE UNIFORMS in " << name() << " " << infologLength << std::endl;
	// glGetProgramiv(id_, GL_INFO_LOG_LENGTH, &infologLength);
	// if (infologLength > 1)
	// {
	// 	char* infoLog = new char[infologLength];
	// 	int charsWritten = 0;
	// 	glGetProgramInfoLog(id_, infologLength, &charsWritten, infoLog);
	// 	std::cerr << "Link message: " << infoLog << std::endl;
	// 	delete[] infoLog;
	// }

	get_matrices_uniforms();
}

/*****************************************************************************/
// ShaderParam
/*****************************************************************************/

ShaderParam::ShaderParam(ShaderProgram* prg) : shader_(prg), attributes_initialized_(false)
{
	vao_ = std::make_unique<VAO>();
}

void ShaderParam::bind(const GLMat4& proj, const GLMat4& mv)
{
	shader_->bind();
	shader_->set_matrices(proj, mv);
	bind_texture_buffers();
	set_uniforms();
	vao_->bind();
}

void ShaderParam::bind()
{
	shader_->bind();
	bind_texture_buffers();
	set_uniforms();
	vao_->bind();
}

void ShaderParam::release()
{
	vao_->release();
	release_texture_buffers();
	shader_->release();
}

void ShaderParam::bind_texture_buffers()
{
}

void ShaderParam::release_texture_buffers()
{
}

void ShaderParam::set_texture_buffer_vbo(uint32, VBO*)
{
}

void ShaderParam::set_vbos(const std::vector<VBO*>& vbos)
{
	assert(uint32(vbos.size()) == shader_->nb_attributes());

	if (shader_->use_texture_buffer())
	{
		attributes_initialized_ = true;
		for (uint32 i = 0; i < uint32(vbos.size()); ++i)
		{
			set_texture_buffer_vbo(i, vbos[i]);
			if (!vbos[i])
				attributes_initialized_ = false;
		}
		return;
	}
	else
	{
		attributes_initialized_ = true;
		vao_->bind();
		GLuint attrib = 1u;
		for (auto* v : vbos)
		{
			if (v)
				v->associate(attrib);
			else
				attributes_initialized_ = false;
			attrib++;
		}
		vao_->release();
	}
}

} // namespace rendering

} // namespace cgogn

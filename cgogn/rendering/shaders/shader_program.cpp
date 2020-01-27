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

#include <cgogn/rendering/shaders/shader_program.h>

namespace cgogn
{

namespace rendering
{

const GLColor ShaderParam::color_front_default = GLColor(0, 0.8f, 0, 1);
const GLColor ShaderParam::color_back_default = GLColor(0, 0, 0.8f, 1);
const GLColor ShaderParam::color_ambiant_default = GLColor(0.1f, 0.1f, 0, 1);
const GLColor ShaderParam::color_spec_default = GLColor(1, 1, 1, 1);
const GLColor ShaderParam::color_line_default = GLColor(1, 1, 0, 1);
const GLColor ShaderParam::color_point_default = GLColor(1, 1, 1, 1);

void Shader::compile(const std::string& src)
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
				  << "compilation de "
				  << "msg"
				  << " : " << std::endl
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
			while ((a < line.size()) && (line[a] >= '0') && (line[a] <= '9'))
				a++;
			std::size_t b = a + 1;
			while ((b < line.size()) && (line[b] >= '0') && (line[b] <= '9'))
				b++;
			if (b < line.size())
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
	GL_ASSERT()
}

ShaderProgram::ShaderProgram() : vert_shader_(nullptr), frag_shader_(nullptr), geom_shader_(nullptr)
{
	id_ = glCreateProgram();  GL_ASSERT()
}

// ShaderProgram::ShaderProgram(const std::string& vert_src, const std::string& frag_src):
//	vert_shader_(nullptr),
//	frag_shader_(nullptr),
//	geom_shader_(nullptr)
//{
//	id_ = glCreateProgram();

//	load(vert_src, frag_src);
//}

// ShaderProgram::ShaderProgram(const std::string& vert_src, const std::string& frag_src, const std::string& geom_src):
//	vert_shader_(nullptr),
//	frag_shader_(nullptr),
//	geom_shader_(nullptr)
//{
//	id_ = glCreateProgram();

//	load(vert_src, frag_src, geom_src);
//}

ShaderProgram::~ShaderProgram()
{
	if (vert_shader_)
		delete vert_shader_;
	if (geom_shader_)
		delete geom_shader_;
	if (frag_shader_)
		delete frag_shader_;

	glDeleteProgram(id_);
	GL_ASSERT()
}

void ShaderProgram::load(const std::string& vert_src, const std::string& frag_src)
{
	vert_shader_ = new Shader(GL_VERTEX_SHADER);
	vert_shader_->compile(vert_src);

	frag_shader_ = new Shader(GL_FRAGMENT_SHADER);
	frag_shader_->compile(frag_src);

	glAttachShader(id_, vert_shader_->shaderId());
	glAttachShader(id_, frag_shader_->shaderId());

	glLinkProgram(id_);

	// puis detache (?)
	glDetachShader(id_, frag_shader_->shaderId());
	glDetachShader(id_, vert_shader_->shaderId());

	// Print log if needed
	int infologLength = 0;
	glGetProgramiv(id_, GL_INFO_LOG_LENGTH, &infologLength);
	if (infologLength > 1)
	{
		char* infoLog = new char[infologLength];
		int charsWritten = 0;
		glGetProgramInfoLog(id_, infologLength, &charsWritten, infoLog);
		std::cerr << "Link message: " << infoLog << std::endl;
		delete[] infoLog;
	}

	get_matrices_uniforms();
	GL_ASSERT()
}

void ShaderProgram::load(const std::string& vert_src, const std::string& frag_src, const std::string& geom_src)
{
	vert_shader_ = new Shader(GL_VERTEX_SHADER);
	vert_shader_->compile(vert_src);

	geom_shader_ = new Shader(GL_GEOMETRY_SHADER);
	geom_shader_->compile(geom_src);

	frag_shader_ = new Shader(GL_FRAGMENT_SHADER);
	frag_shader_->compile(frag_src);

	glAttachShader(id_, vert_shader_->shaderId());
	glAttachShader(id_, geom_shader_->shaderId());
	glAttachShader(id_, frag_shader_->shaderId());

	glLinkProgram(id_);

	// puis detache (?)
	glDetachShader(id_, frag_shader_->shaderId());
	glDetachShader(id_, geom_shader_->shaderId());
	glDetachShader(id_, vert_shader_->shaderId());

	// Print log if needed
	GLint infologLength = 0;
	glGetProgramiv(id_, GL_LINK_STATUS, &infologLength);
	if (infologLength != GL_TRUE)
		std::cerr << "PB GL_LINK_STATUS" << std::endl;
	glGetProgramiv(id_, GL_VALIDATE_STATUS, &infologLength);
	if (infologLength != GL_TRUE)
		std::cerr << "PB GL_VALIDATE_STATUS" << std::endl;

	glGetProgramiv(id_, GL_INFO_LOG_LENGTH, &infologLength);
	if (infologLength > 1)
	{
		char* infoLog = new char[infologLength];
		int charsWritten = 0;
		glGetProgramInfoLog(id_, infologLength, &charsWritten, infoLog);
		std::cerr << "Link message: " << infoLog << std::endl;
		delete[] infoLog;
	}

	get_matrices_uniforms();
	GL_ASSERT()
}

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

void ShaderProgram::get_matrices_uniforms()
{
	unif_mvp_matrix_ = -1;
	unif_mv_matrix_ = -1;
	unif_projection_matrix_ = -1;
	unif_normal_matrix_ = -1;
	unif_mvp_matrix_ = glGetUniformLocation(id_, "mvp_matrix");
	unif_mv_matrix_ = glGetUniformLocation(id_, "model_view_matrix");
	unif_projection_matrix_ = glGetUniformLocation(id_, "projection_matrix");
	unif_normal_matrix_ = glGetUniformLocation(id_, "normal_matrix");  GL_ASSERT()
}

void ShaderProgram::set_matrices(const GLMat4d& proj, const GLMat4d& mv)
{
	if (unif_mvp_matrix_ >= 0)
	{
		GLMat4d mvp = proj * mv;
		GLMat4 m = mvp.cast<float>();
		glUniformMatrix4fv(unif_mvp_matrix_, 1, false, m.data());  GL_ASSERT()
	}
	if (unif_projection_matrix_ >= 0)
	{
		GLMat4 m = proj.cast<float>();
		glUniformMatrix4fv(unif_projection_matrix_, 1, false, m.data());  GL_ASSERT()
	}
	if (unif_mv_matrix_ >= 0)
	{
		GLMat4 m = mv.cast<float>();
		glUniformMatrix4fv(unif_mv_matrix_, 1, false, m.data());  GL_ASSERT()
	}
	if (unif_normal_matrix_ >= 0)
	{
		Eigen::Affine3d t(mv);
		GLMat3 normal_matrix = t.linear().inverse().transpose().matrix().cast<float>();
		glUniformMatrix3fv(unif_normal_matrix_, 1, false, normal_matrix.data());  GL_ASSERT()
	}
}

void ShaderProgram::set_matrices(const GLMat4& proj, const GLMat4& mv)
{
	if (unif_mvp_matrix_ >= 0)
	{
		GLMat4 m = proj * mv;
		glUniformMatrix4fv(unif_mvp_matrix_, 1, false, m.data());  GL_ASSERT()
	}
	if (unif_projection_matrix_ >= 0)
		glUniformMatrix4fv(unif_projection_matrix_, 1, false, proj.data());  GL_ASSERT()
	if (unif_mv_matrix_ >= 0)
		glUniformMatrix4fv(unif_mv_matrix_, 1, false, mv.data());  GL_ASSERT()
	if (unif_normal_matrix_ >= 0)
	{
		Eigen::Affine3d t(mv.cast<float64>());
		GLMat3 normal_matrix = t.linear().inverse().transpose().matrix().cast<float32>();
		glUniformMatrix3fv(unif_normal_matrix_, 1, false, normal_matrix.data());  GL_ASSERT()
	}
}

void ShaderProgram::set_view_matrix(const GLMat4d& mv)
{
	if (unif_mv_matrix_ >= 0)
	{
		GLMat4 m = mv.cast<float>();
		glUniformMatrix4fv(unif_mv_matrix_, 1, false, m.data());  GL_ASSERT()
	}
	if (unif_normal_matrix_ >= 0)
	{
		Eigen::Affine3d t(mv);
		GLMat3 normal_matrix = t.linear().inverse().transpose().matrix().cast<float32>();
		glUniformMatrix3fv(unif_normal_matrix_, 1, false, normal_matrix.data());  GL_ASSERT()
	}
}

void ShaderProgram::set_view_matrix(const GLMat4& mv)
{
	if (unif_mv_matrix_ >= 0)
		glUniformMatrix4fv(unif_mv_matrix_, 1, false, mv.data());  GL_ASSERT()
	if (unif_normal_matrix_ >= 0)
	{
		Eigen::Affine3d t(mv.cast<float64>());
		GLMat3 normal_matrix = t.linear().inverse().transpose().matrix().cast<float>();
		glUniformMatrix3fv(unif_normal_matrix_, 1, false, normal_matrix.data());  GL_ASSERT()
	}
}

ShaderParam::ShaderParam(ShaderProgram* prg) : shader_(prg), vao_initialized_(false)
{
	vao_ = std::make_unique<VAO>();
	vao_->create();  GL_ASSERT()
	vao_initialized_ = true;
}

void ShaderParam::bind(const GLMat4& proj, const GLMat4& mv)
{
	shader_->bind();  GL_ASSERT()
	shader_->set_matrices(proj, mv);  GL_ASSERT()
	set_uniforms();  GL_ASSERT()
	vao_->bind();  GL_ASSERT()
}

void ShaderParam::bind()
{
	shader_->bind();  GL_ASSERT()
	set_uniforms();  GL_ASSERT()
	vao_->bind();  GL_ASSERT()
}

void ShaderParam::release()
{
	vao_->release();  GL_ASSERT()
	shader_->release();  GL_ASSERT()
}

} // namespace rendering

} // namespace cgogn

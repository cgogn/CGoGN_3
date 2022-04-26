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

#ifndef CGOGN_RENDERING_SHADERS_SHADER_PROGRAM_H_
#define CGOGN_RENDERING_SHADERS_SHADER_PROGRAM_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/rendering/types.h>
#include <cgogn/rendering/vao.h>

#include <array>
#include <iostream>
#include <memory>

#define DECLARE_SHADER_CLASS(NAME, TB, STRNAME)                                                                        \
	class ShaderParam##NAME;                                                                                           \
	class CGOGN_RENDERING_EXPORT Shader##NAME : public ShaderProgram                                                   \
	{                                                                                                                  \
	public:                                                                                                            \
		using Self = Shader##NAME;                                                                                     \
		using Param = ShaderParam##NAME;                                                                               \
		friend Param;                                                                                                  \
		inline static std::unique_ptr<Param> generate_param()                                                          \
		{                                                                                                              \
			if (!instance_)                                                                                            \
			{                                                                                                          \
				instance_ = new Self();                                                                                \
				ShaderProgram::register_instance(instance_);                                                           \
			}                                                                                                          \
			return std::make_unique<Param>(instance_);                                                                 \
		}                                                                                                              \
		inline std::string name() const override                                                                       \
		{                                                                                                              \
			return STRNAME;                                                                                            \
		}                                                                                                              \
		inline bool use_texture_buffer() const override                                                                \
		{                                                                                                              \
			return TB;                                                                                                 \
		}                                                                                                              \
                                                                                                                       \
	protected:                                                                                                         \
		Shader##NAME();                                                                                                \
		Shader##NAME(const Shader##NAME&) = delete;                                                                    \
		Shader##NAME(Shader##NAME&&) = delete;                                                                         \
		Shader##NAME& operator=(const Shader##NAME&) = delete;                                                         \
		Shader##NAME& operator=(Shader##NAME&&) = delete;                                                              \
		static Self* instance_;                                                                                        \
	};

namespace cgogn
{

namespace rendering
{

inline GLColor color(uint8 R, uint8 G, uint8 B, uint8 A = 255u)
{
	return GLColor(float32(R) / 255.0f, float32(G) / 255.0f, float32(B) / 255.0f, float32(A) / 255.0f);
}

// convenient conversion function
inline void* void_ptr(uint32 x)
{
	return reinterpret_cast<void*>(uint64_t(x));
}

/*****************************************************************************/
// Shader
/*****************************************************************************/

class CGOGN_RENDERING_EXPORT Shader
{
protected:
	GLuint id_;

public:
	Shader() = delete;
	Shader(const Shader&) = delete;
	Shader& operator=(const Shader&) = delete;

	inline Shader(GLenum type)
	{
		id_ = glCreateShader(type);
	}

	inline ~Shader()
	{
		glDeleteShader(id_);
	}

	inline GLuint id() const
	{
		return id_;
	}

	void compile(const std::string& src, const std::string& prg_name);
};

/*****************************************************************************/
// ShaderProgram
/*****************************************************************************/

class CGOGN_RENDERING_EXPORT ShaderProgram
{
protected:
	static std::vector<ShaderProgram*>* instances_;

	GLuint id_;
	Shader* vertex_shader_;
	Shader* fragment_shader_;
	Shader* geometry_shader_;

	GLint uniform_mvp_matrix_;
	GLint uniform_mv_matrix_;
	GLint uniform_projection_matrix_;
	GLint uniform_normal_matrix_;

	uint32 nb_attributes_;

	std::vector<GLint> uniforms_;

public:
	ShaderProgram();
	ShaderProgram(const ShaderProgram&) = delete;
	ShaderProgram& operator=(const ShaderProgram&) = delete;

	virtual ~ShaderProgram();

	static void register_instance(ShaderProgram* sh);
	static void clean_all();

	inline GLuint id() const
	{
		return id_;
	}

	virtual std::string name() const = 0;

	virtual bool use_texture_buffer() const = 0;

	inline uint32 nb_attributes() const
	{
		return nb_attributes_;
	}

	inline void bind()
	{
		glUseProgram(id_);
	}

	inline void release()
	{
		glUseProgram(0);
	}

	inline void get_uniform(const GLchar* str)
	{
		GLint u = glGetUniformLocation(id_, str);
		if (u >= 0)
			uniforms_.push_back(u);
		else
			std::cerr << "Warning uniform " << str << " does not exist in shader " << name() << std::endl;
	}

	template <typename T1>
	void get_uniforms(T1 p1)
	{
		get_uniform(p1);
	}

	template <typename T1, typename... Ts>
	void get_uniforms(T1 p1, Ts... pn)
	{
		get_uniform(p1);
		get_uniforms(pn...);
	}

	inline void set_uniform_value(std::size_t i, const float32 v)
	{
		glUniform1f(uniforms_[i], v);
	}
	inline void set_uniform_value(std::size_t i, const GLVec2& v)
	{
		glUniform2fv(uniforms_[i], 1, v.data());
	}
	inline void set_uniform_value(std::size_t i, const GLVec3& v)
	{
		glUniform3fv(uniforms_[i], 1, v.data());
	}
	inline void set_uniform_value(std::size_t i, const GLVec4& v)
	{
		glUniform4fv(uniforms_[i], 1, v.data());
	}
	inline void set_uniform_value(std::size_t i, const int32 v)
	{
		glUniform1i(uniforms_[i], v);
	}
	inline void set_uniform_value(std::size_t i, const uint32 v)
	{
		glUniform1ui(uniforms_[i], v);
	}
	inline void set_uniform_value(std::size_t i, const bool v)
	{
		glUniform1i(uniforms_[i], int32(v));
	}

	template <typename T>
	void set_uniforms_values(T v)
	{
		set_uniform_value(uint32(uniforms_.size()) - 1, v);
	}

	template <typename T, typename... Ts>
	void set_uniforms_values(T v, Ts... vs)
	{
		set_uniform_value(uint32(uniforms_.size()) - 1 - sizeof...(Ts), v);
		set_uniforms_values(vs...);
	}

	void get_matrices_uniforms();

	void set_matrices(const GLMat4& proj, const GLMat4& mv);
	void set_matrices(const GLMat4d& proj, const GLMat4d& mv);

	inline void bind_attrib_location(GLuint attrib, const char* str_var)
	{
		glBindAttribLocation(id_, attrib, str_var);
	}

	inline void bind_attrib_location(GLuint attrib, const std::string& str_var)
	{
		glBindAttribLocation(id_, attrib, str_var.c_str());
	}

	template <typename T1>
	void internal_bind_attrib_locations(GLuint attrib1, T1 p1)
	{
		bind_attrib_location(attrib1, p1);
	}

	template <typename T1, typename... Ts>
	void internal_bind_attrib_locations(GLuint attrib1, T1 p1, Ts... pn)
	{
		bind_attrib_location(attrib1, p1);
		internal_bind_attrib_locations(attrib1 + 1, pn...);
	}

	void bind_attrib_locations()
	{
	}

	template <typename... Ts>
	void bind_attrib_locations(Ts... pn)
	{
		internal_bind_attrib_locations(1u, pn...);
	}

	void load(const std::string& vertex_program_src, const std::string& fragment_program_src);
	void load3(const std::string& vertex_program_src, const std::string& fragment_program_src,
			   const std::string& geometry_program_src);

	template <typename... Ts>
	void load3_bind(const std::string& vertex_program_src, const std::string& fragment_program_src,
					const std::string& geometry_program_src, Ts... pn)
	{
		vertex_shader_ = new Shader(GL_VERTEX_SHADER);
		vertex_shader_->compile(vertex_program_src, name());

		geometry_shader_ = new Shader(GL_GEOMETRY_SHADER);
		geometry_shader_->compile(geometry_program_src, name());

		fragment_shader_ = new Shader(GL_FRAGMENT_SHADER);
		fragment_shader_->compile(fragment_program_src, name());

		glAttachShader(id_, vertex_shader_->id());
		glAttachShader(id_, geometry_shader_->id());
		glAttachShader(id_, fragment_shader_->id());

		nb_attributes_ = sizeof...(Ts);
		bind_attrib_locations(pn...);

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
		// 	std::cerr << "PB GL_LINK_STATUS load3_bind " << name() << std::endl;
		// glGetProgramiv(id_, GL_VALIDATE_STATUS, &infologLength);
		// if (infologLength != GL_TRUE)
		// 	std::cerr << "PB GL_VALIDATE_STATUS load3_bind " << name() << std::endl;
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

	template <typename... Ts>
	void load2_bind(const std::string& vertex_program_src, const std::string& fragment_program_src, Ts... pn)
	{
		vertex_shader_ = new Shader(GL_VERTEX_SHADER);
		vertex_shader_->compile(vertex_program_src, name());

		fragment_shader_ = new Shader(GL_FRAGMENT_SHADER);
		fragment_shader_->compile(fragment_program_src, name());

		glAttachShader(id_, vertex_shader_->id());
		glAttachShader(id_, fragment_shader_->id());

		nb_attributes_ = sizeof...(Ts);
		bind_attrib_locations(pn...);

		glLinkProgram(id_);

		// puis detache (?)
		glDetachShader(id_, fragment_shader_->id());
		glDetachShader(id_, vertex_shader_->id());

		// glValidateProgram(id_);
		// // Print log if needed
		// GLint infologLength = 0;
		// glGetProgramiv(id_, GL_LINK_STATUS, &infologLength);
		// if (infologLength != GL_TRUE)
		// 	std::cerr << "PB GL_LINK_STATUS load2_bind " << name() << " " << infologLength << std::endl;
		// glGetProgramiv(id_, GL_VALIDATE_STATUS, &infologLength);
		// if (infologLength != GL_TRUE)
		// 	std::cerr << "PB GL_VALIDATE_STATUS load2_bind " << name() << " " << infologLength << std::endl;
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

	template <typename... Ts>
	void load_tfb1_bind(const std::string& vertex_program_src, const std::vector<std::string>& tf_outs, Ts... pn)
	{
		vertex_shader_ = new Shader(GL_VERTEX_SHADER);
		vertex_shader_->compile(vertex_program_src, name());

		glAttachShader(id_, vertex_shader_->id());

		nb_attributes_ = sizeof...(Ts);
		bind_attrib_locations(pn...);

		if (!tf_outs.empty())
		{
			std::vector<const char*> tfo;
			for (const auto& t : tf_outs)
				tfo.push_back(t.c_str());
			glTransformFeedbackVaryings(id_, GLsizei(uint32(tf_outs.size())), tfo.data(), GL_SEPARATE_ATTRIBS);
		}

		glLinkProgram(id_);

		glDetachShader(id_, vertex_shader_->id());

		// glValidateProgram(id_);
		// // Print log if needed
		// GLint infologLength = 0;
		// glGetProgramiv(id_, GL_LINK_STATUS, &infologLength);
		// if (infologLength != GL_TRUE)
		// 	std::cerr << "PB GL_LINK_STATUS load_tfb1_bind " << name() << std::endl;
		// glGetProgramiv(id_, GL_VALIDATE_STATUS, &infologLength);
		// if (infologLength != GL_TRUE)
		// 	std::cerr << "PB GL_VALIDATE_STATUS load_tfb1_bind " << name() << std::endl;
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
};

/*****************************************************************************/
// ShaderParam
/*****************************************************************************/

class CGOGN_RENDERING_EXPORT ShaderParam
{
protected:
	ShaderProgram* shader_;
	std::unique_ptr<VAO> vao_;
	bool attributes_initialized_;
	bool optional_clipping_attribute_;

	virtual void set_uniforms() = 0;

	virtual void bind_texture_buffers();
	virtual void release_texture_buffers();
	virtual void set_texture_buffer_vbo(uint32, VBO*);

public:
	ShaderParam(ShaderProgram* prg, bool opt_clip = false);
	ShaderParam(const ShaderParam&) = delete;
	ShaderParam& operator=(const ShaderParam&) = delete;

	inline virtual ~ShaderParam()
	{
	}

	inline bool attributes_initialized() const
	{
		return attributes_initialized_; // >= (1u << shader_->nb_attributes()) - 1;
	}

	inline void set_vao_name(const std::string& name)
	{
		vao_->set_name(name);
	}

	/**
	 * @brief bind the shader, set uniforms & matrices, bind vao
	 * @param proj projection matrix
	 * @param mv modelview matrix
	 */
	void bind(const GLMat4& proj, const GLMat4& mv);

	void bind();

	/**
	 * @brief release vao and shader
	 */
	void release();

	/**
	 * @brief set vbos into the vao
	 * @param all vbos in order of attribs
	 */
	virtual void set_vbos(const std::vector<VBO*>& vbos);

	// /**
	//  * @brief set one vbo into the vao
	//  * @param attrib_id, vbo
	//  */
	// virtual void set_vbo(GLuint att, VBO* vbo);
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_SHADER_PROGRAM_H_

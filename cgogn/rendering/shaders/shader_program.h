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

#ifndef CGOGN_RENDERING_SHADERS_SHADERPROGRAM_H_
#define CGOGN_RENDERING_SHADERS_SHADERPROGRAM_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/types.h>
#include <cgogn/rendering/vao.h>

#include <iostream>
#include <memory>

#define DECLARE_SHADER_CLASS(NAME, STRNA)                                                                              \
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
			return STRNA;                                                                                              \
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

struct PossibleParameters
{
	GLColor color_;
	GLColor ambiant_color_;
	GLColor front_color_;
	GLColor back_color_;
	GLColor specular_color_;
	float32 specular_coef_;
	GLVec3 light_position_;
	bool double_side_;
	float32 width_;
	float32 size_;
	float32 explode_;
	bool blending_;
	GLVec4 plane_clip_;
	GLVec4 plane_clip2_;
};

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

	inline GLuint shaderId() const
	{
		return id_;
	}

	void compile(const std::string& src, const std::string& prg_name);
};

class CGOGN_RENDERING_EXPORT ShaderProgram
{
protected:
	static std::vector<ShaderProgram*>* instances_;
	GLuint id_;
	Shader* vert_shader_;
	Shader* frag_shader_;
	Shader* geom_shader_;

	void load(const std::string& vert_src, const std::string& frag_src);
	void load(const std::string& vert_src, const std::string& frag_src, const std::string& geom_src);

	GLint unif_mvp_matrix_;
	GLint unif_mv_matrix_;
	GLint unif_projection_matrix_;
	GLint unif_normal_matrix_;
	uint32 nb_attributes_;

	std::vector<GLint> uniforms_;

public:
	//	enum Attrib_Indices: GLuint
	//	{
	//		ATTRIB_POS     = 1u,
	//		ATTRIB_COLOR   = 2u,
	//		ATTRIB_NORM    = 3u,
	//		ATTRIB_TC      = 4u,
	//		ATTRIB_SIZE    = 5u,
	//		ATTRIB_CUSTOM1 = 5u,
	//		ATTRIB_CUSTOM2 = 6u,
	//		ATTRIB_CUSTOM3 = 7u,
	//		ATTRIB_CUSTOM4 = 8u,
	//	};

	static void register_instance(ShaderProgram* sh);

	static void clean_all();

	ShaderProgram();
	ShaderProgram(const ShaderProgram&) = delete;
	ShaderProgram& operator=(const ShaderProgram&) = delete;

	virtual ~ShaderProgram();

	inline GLuint id() const
	{
		return id_;
	}

	//	inline void start_use()
	//	{
	//		glUseProgram(id_);
	//	}

	//	inline void stop_use()
	//	{
	//		glUseProgram(0);
	//	}

	virtual std::string name() const = 0;

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

	inline GLint uniform_location(const GLchar* str) const
	{
		return glGetUniformLocation(id_, str);
	}

	inline void add_uniform(const GLchar* str)
	{
		GLint u = glGetUniformLocation(id_, str);
		if (u >= 0)
			uniforms_.push_back(u);
		else
			std::cerr << "Warning uniform " << str << " does not exist in shader " << name() << std::endl;
	}

	template <typename T1>
	void add_uniforms(T1 p1)
	{
		add_uniform(p1);
	}

	template <typename T1, typename... Ts>
	void add_uniforms(T1 p1, Ts... pn)
	{
		add_uniform(p1);
		add_uniforms(pn...);
	}

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

	template <typename T1>
	void set_uniforms_values(T1 p1)
	{
		set_uniform_value(uniforms_.size() - 1, p1);
	}

	template <typename T1, typename... Ts>
	void set_uniforms_values(T1 p1, Ts... pn)
	{
		set_uniform_value(uniforms_.size() - 1 - sizeof...(Ts), p1);
		set_uniforms_values(pn...);
	}

	void get_matrices_uniforms();

	void set_matrices(const GLMat4& proj, const GLMat4& mv);
	void set_matrices(const GLMat4d& proj, const GLMat4d& mv);

	void set_view_matrix(const GLMat4& mv);
	void set_view_matrix(const GLMat4d& mv);

	template <typename... Ts>
	void load3_bind(const std::string& vert_src, const std::string& frag_src, const std::string& geom_src, Ts... pn)
	{
		vert_shader_ = new Shader(GL_VERTEX_SHADER);
		vert_shader_->compile(vert_src, name());

		geom_shader_ = new Shader(GL_GEOMETRY_SHADER);
		geom_shader_->compile(geom_src, name());

		frag_shader_ = new Shader(GL_FRAGMENT_SHADER);
		frag_shader_->compile(frag_src, name());

		glAttachShader(id_, vert_shader_->shaderId());
		glAttachShader(id_, geom_shader_->shaderId());
		glAttachShader(id_, frag_shader_->shaderId());

		//		set_locations();
		nb_attributes_ = sizeof...(Ts);
		bind_attrib_locations(pn...);

		glLinkProgram(id_);

		// puis detache (?)
		glDetachShader(id_, frag_shader_->shaderId());
		glDetachShader(id_, geom_shader_->shaderId());
		glDetachShader(id_, vert_shader_->shaderId());

		glValidateProgram(id_);
		// Print log if needed
		GLint infologLength = 0;
		glGetProgramiv(id_, GL_LINK_STATUS, &infologLength);
		if (infologLength != GL_TRUE)
			std::cerr << "PB GL_LINK_STATUS load3_bind " << name() << std::endl;
		glGetProgramiv(id_, GL_VALIDATE_STATUS, &infologLength);
		if (infologLength != GL_TRUE)
			std::cerr << "PB GL_VALIDATE_STATUS load3_bind " << name() << std::endl;

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
	}

	template <typename... Ts>
	void load2_bind(const std::string& vert_src, const std::string& frag_src, Ts... pn)
	{
		vert_shader_ = new Shader(GL_VERTEX_SHADER);
		vert_shader_->compile(vert_src, name());

		frag_shader_ = new Shader(GL_FRAGMENT_SHADER);
		frag_shader_->compile(frag_src, name());

		glAttachShader(id_, vert_shader_->shaderId());
		glAttachShader(id_, frag_shader_->shaderId());

		nb_attributes_ = sizeof...(Ts);
		bind_attrib_locations(pn...);

		glLinkProgram(id_);

		// puis detache (?)
		glDetachShader(id_, frag_shader_->shaderId());
		glDetachShader(id_, vert_shader_->shaderId());

		glValidateProgram(id_);
		// Print log if needed
		GLint infologLength = 0;
		glGetProgramiv(id_, GL_LINK_STATUS, &infologLength);
		if (infologLength != GL_TRUE)
			std::cerr << "PB GL_LINK_STATUS load2_bind " << name() << " " << infologLength << std::endl;
		glGetProgramiv(id_, GL_VALIDATE_STATUS, &infologLength);
		if (infologLength != GL_TRUE)
			std::cerr << "PB GL_VALIDATE_STATUS load2_bind " << name() << " " << infologLength << std::endl;

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
	}

	template <typename... Ts>
	void load_tfb1_bind(const std::string& vert_src, const std::vector<std::string>& tf_outs, Ts... pn)
	{
		vert_shader_ = new Shader(GL_VERTEX_SHADER);
		vert_shader_->compile(vert_src, name());

		glAttachShader(id_, vert_shader_->shaderId());

		nb_attributes_ = sizeof...(Ts);
		bind_attrib_locations(pn...);

		if (!tf_outs.empty())
		{
			std::vector<const char*> tfo;
			for (const auto& t : tf_outs)
				tfo.push_back(t.c_str());
			glTransformFeedbackVaryings(id_, GLsizei(tf_outs.size()), tfo.data(), GL_SEPARATE_ATTRIBS);
		}

		glLinkProgram(id_);

		glDetachShader(id_, vert_shader_->shaderId());

		glValidateProgram(id_);
		// Print log if needed
		GLint infologLength = 0;
		glGetProgramiv(id_, GL_LINK_STATUS, &infologLength);
		if (infologLength != GL_TRUE)
			std::cerr << "PB GL_LINK_STATUS load_tfb1_bind " << name() << std::endl;
		glGetProgramiv(id_, GL_VALIDATE_STATUS, &infologLength);
		if (infologLength != GL_TRUE)
			std::cerr << "PB GL_VALIDATE_STATUS load_tfb1_bind " << name() << std::endl;
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
	}
};

class CGOGN_RENDERING_EXPORT ShaderParam
{
protected:
	ShaderProgram* shader_;
	std::unique_ptr<VAO> vao_;
	bool vao_initialized_;

	virtual void set_uniforms() = 0;

	static const GLColor color_front_default;
	static const GLColor color_back_default;
	static const GLColor color_ambiant_default;
	static const GLColor color_spec_default;
	static const GLColor color_line_default;
	static const GLColor color_point_default;

public:
	ShaderParam(ShaderProgram* prg);
	ShaderParam(const ShaderParam&) = delete;
	ShaderParam& operator=(const ShaderParam&) = delete;

	virtual void pick_parameters(const PossibleParameters&);

	inline virtual ~ShaderParam()
	{
	}

	inline bool vao_initialized() const
	{
		return vao_initialized_;
	}

	inline void bind_vao()
	{
		vao_->bind();
	}

	inline void release_vao()
	{
		vao_->release();
	}

	inline ShaderProgram* get_shader()
	{
		return shader_;
	}

	inline void set_name(const std::string& name)
	{
		vao_->set_name(name);
	}

	/**
	 * @brief bind the shader set uniforms & matrices, bind vao
	 * @param proj projection matrix
	 * @param mv model-view matrix
	 */
	void bind(const GLMat4& proj, const GLMat4& mv);

	void bind();

	/**
	 * @brief release vao and shader
	 */
	void release();

	/**
	 * @brief set vbos into the vao
	 * @param vbos
	 */
	virtual void set_vbos(const std::vector<VBO*>& vbos);
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHADERS_SHADERPROGRAM_H_

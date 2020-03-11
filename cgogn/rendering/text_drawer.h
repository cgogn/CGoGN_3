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

#ifndef CGOGN_RENDERING_TEXT_DRAWER_H_
#define CGOGN_RENDERING_TEXT_DRAWER_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shaders/shader_text.h>

#include <array>

namespace cgogn
{

namespace rendering
{

// class CGOGN_RENDERING_EXPORT TextDrawerEnd {};

/**
 * @brief Rendering of volumes
 *
 * Typical usage:
 *
 *  std::unique_ptr<cgogn::rendering::TextDrawer> text_;	// can be shared between contexts
 *  std::unique_ptr<cgogn::rendering::TextDrawer::Renderer> text_rend_; // one by context,
 *
 * init:
 *  text_ = cgogn::make_unique<cgogn::rendering::TextDrawer>();
 *  text_rend_ = text_->generate_renderer();
 *
 *  *text_ << pos << size << color << str; // add str at pos using size an color
 *  *text_ << pos << color << str;         // add str at pos using color (keep size)
 *  *text_ << pos << size << str;          // add str at pos using size (keep color)
 *  *text_ << pos << str;                  // add str at pos (keep size & color)
 *  *text_ << str;                         // concat str with previous string
 *  text_->begin();
 * draw:
 *  text_rend_->draw(proj, view, this);
 *
 */
class CGOGN_RENDERING_EXPORT TextDrawer
{
protected:
	using Vec3f = std::array<float32, 3>;
	using Vec4f = std::array<float32, 4>;

	std::unique_ptr<VBO> vbo_pos_;
	std::unique_ptr<VBO> vbo_char_;
	std::unique_ptr<VBO> vbo_colsz_;

	static Texture2D* texture_;

	std::vector<Vec3f> positions_;
	std::vector<std::string> strings_;
	std::vector<GLColor> colors_;
	std::vector<float32> sizes_;

	Vec3f current_pos_;
	GLColor current_color_;
	float32 current_size_;
	bool next_pos_;

public:
	using Self = TextDrawer;

	/**
	 * constructor, init all buffers (data and OpenGL) and shader
	 * @Warning need OpenGL context
	 */
	TextDrawer();

	/**
	 * release buffers and shader
	 */
	~TextDrawer();

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(TextDrawer);

	class End
	{
	};

	static TextDrawer::End end;

	Self& operator<<(TextDrawer::End);

	template <typename VEC3>
	inline auto operator<<(const VEC3& p) -> typename std::enable_if<geometry::vector_traits<VEC3>::OK, Self&>::type
	{
		current_pos_ = Vec3f{float32(p[0]), float32(p[1]), float32(p[2])};
		next_pos_ = true;
		return *this;
	}

	Self& operator<<(const GLColor& col);

	Self& operator<<(float32 sz);

	Self& operator<<(const std::string& str);

	inline Self& operator<<(const char* cstr)
	{
		return this->operator<<(std::string(cstr));
	}

	void update_text(std::size_t pos, const std::string& str);

	void scale_text(float sc);

	class CGOGN_RENDERING_EXPORT Renderer
	{
		friend class TextDrawer;
		std::unique_ptr<ShaderText::Param> param_text_;
		TextDrawer* text_drawer_data_;
		Renderer(TextDrawer* td);

	public:
		~Renderer();

		void draw(const GLMat4& projection, const GLMat4& modelview);

		/**
		 * @brief set italic %
		 * @param it [0, 0.5]
		 */
		void inline set_italic(float32 it)
		{
			param_text_->italic_ = it;
		}
	};

	/**
	 * @brief generate a renderer (one per context)
	 * @return pointer on renderer
	 */
	inline std::unique_ptr<Renderer> generate_renderer()
	{
		return std::unique_ptr<Renderer>(new Renderer(this));
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_TEXT_DRAWER_H_

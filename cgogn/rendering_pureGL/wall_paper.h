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

#ifndef CGOGN_RENDERING_WALL_PAPER_H_
#define CGOGN_RENDERING_WALL_PAPER_H_

#include <cgogn/rendering_pureGL/shaders/shader_texture.h>
#include <cgogn/rendering_pureGL/vbo.h>
#include <cgogn/rendering_pureGL/cgogn_rendering_puregl_export.h>

namespace cgogn
{

namespace rendering_pgl
{

/**
 * @brief WallPaper: allow rendering of a texture, front or back, full-screen or not
 *
 * Typical usage:
 *
 *  std::unique_ptr<cgogn::rendering::WallPaper> wp_;	// can be shared between contexts
 *  std::unique_ptr<cgogn::rendering::WallPaper::Renderer> wp_rend_; // one by context,
 *
 * init:
 *  wp_ = cgogn::make_unique<cgogn::rendering::WallPaper>();
 *  wp_rend_ = wp_->generate_renderer();
 *  wp_->update<Vec3>(map_,vertex_position_);
 *
 * draw:
 *  wp_rend_->draw(proj,view,this);
 *
 */
class CGOGN_RENDERING_PUREGL_EXPORT WallPaper
{
protected:

	std::unique_ptr<VBO> vbo_pos_;
	std::unique_ptr<VBO> vbo_tc_;
	std::unique_ptr<Texture2D> texture_;
	void init(const GLImage& img);

public:

	class CGOGN_RENDERING_PUREGL_EXPORT Renderer
	{
		friend class WallPaper;

		std::unique_ptr<ShaderTexture::Param> param_texture_;
		WallPaper* wall_paper_data_;

		Renderer(WallPaper* wp);

	public:

		~Renderer();

		void draw();
	};

	using Self = WallPaper;

	/**
	 * @brief constructor, init all buffers (data and OpenGL) and shader
	 * @param img image
	 * @warning need OpenGL context
	 */
	WallPaper(const GLImage& img);

	/**
	 * @brief constructor, init all buffers (data and OpenGL) and shader
	 * @param col color unique color of wallPaper
	 */
	WallPaper(const GLColor& col);

	/**
	 * @brief  constructor, init all buffers (data and OpenGL) and shader
	 * @param col_tl top left color
	 * @param col_tr top right color
	 * @param col_bl bottom left color
	 * @param col_br botton right color
	 */
	WallPaper(const GLColor& col_tl, const GLColor& col_tr, const GLColor& col_bl, const GLColor& col_br);

	/**
	 * @brief change color for unique color image
	 * @param col color
	 */
//	void change_color(const GLColor& col);

	/**
	 * @brief change colors for 4 colors image only
	 * @param col_tl top left color
	 * @param col_tr top right color
	 * @param col_bl bottom left color
	 * @param col_br botton right color
	 */
//	void change_colors(const GLColor& col_tl, const GLColor& col_tr, const GLColor& col_bl, const GLColor& col_br);

	/**
	 * release buffers and shader
	 */
	~WallPaper();

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(WallPaper);

	/**
	 * @brief generate a renderer (one per context)
	 * @return pointer on renderer
	 */
	inline std::unique_ptr<Renderer> generate_renderer()
	{
		return std::unique_ptr<Renderer>(new Renderer(this));
	}

	/**
	 * @brief set the texture in full screen
	 * @param front if true draw with depth of 0 (front) else with depth of ~1 (back)
	 */
	void set_full_screen(bool front = false);

	/**
	 * @brief set a local position for the image in pixel
	 * @warning position & size are converted in % when set, and used even when window size has changed
	 * @param win_w width of window
	 * @param win_h height of window
	 * @param x x pos (0 is left)
	 * @param y y pos (0 is top)
	 * @param w width to draw
	 * @param h height to draw
	 * @param front (default is back drawing)
	 */
	void set_local_position(uint32 win_w, uint32 win_h, uint32 x, uint32 y, uint32 w, uint32 h, bool front = true);

	/**
	 * @brief set a local position for the image in ratio ([0-1]) of the viewport
	 * @param x x pos (0.0 is left)
	 * @param y y pos (0 is top)
	 * @param w width to draw
	 * @param h height to draw
	 * @param front (default is front drawing)
	 */
	void set_local_position(float x, float y, float w, float h, bool front = true);

	void draw();
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_WALL_PAPER_H_

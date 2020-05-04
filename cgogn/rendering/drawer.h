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

#ifndef CGOGN_RENDERING_DRAWER_H_
#define CGOGN_RENDERING_DRAWER_H_

#include <cgogn/rendering/cgogn_rendering_export.h>

#include <cgogn/rendering/shaders/shader_bold_line_color.h>
#include <cgogn/rendering/shaders/shader_color_per_vertex.h>
#include <cgogn/rendering/shaders/shader_flat.h>
#include <cgogn/rendering/shaders/shader_point_sprite.h>
#include <cgogn/rendering/shaders/shader_round_point_color.h>

#include <cgogn/rendering/types.h>

#include <array>

namespace cgogn
{

namespace rendering
{
/**
 * @brief DisplayListDrawer revival of old GL display-list
 *
 * Typical usage:
 *
 *  std::unique_ptr<cgogn::rendering::DisplayListDrawer> drawer_;	// can be shared between contexts
 *  std::unique_ptr<cgogn::rendering::DisplayListDrawer::Renderer> drawer_rend_; // one by context,
 *
 * init:
 *  drawer_ = std::make_unique<cgogn::rendering::DisplayListDrawer>();
 *  drawer_rend_ = drawer_->generate_renderer(); // don't worry automatically deleted when finished
 *  drawer_->new_list();
 *  drawer_->line_width(2.0);
 *  drawer_->begin(GL_LINE_LOOP); // or GL_POINTS, GL_LINES, GL_TRIANGLES
 *    drawer_->color3f(1.0,0.0,0.0);
 *    drawer_->vertex3f(0,0,0);
 *    drawer_->color3f(0.0,1.0,1.0);
 *    ....
 *  drawer_->end();
 *
 * draw:
 *  drawer_rend_->draw(proj,view,this);
 */
class CGOGN_RENDERING_EXPORT DisplayListDrawer
{
	struct PrimParam
	{
		GLint begin;
		GLenum mode;
		float32 width;
		GLsizei nb;
		bool aa;

		PrimParam(std::size_t b, GLenum m, float32 w, bool a) : begin(GLint(b)), mode(m), width(w), nb(0), aa(a)
		{
		}
	};

protected:
	std::unique_ptr<VBO> vbo_pos_;
	std::unique_ptr<VBO> vbo_col_;

	// temporary (between begin()/end()) data storage
	std::vector<GLVec3> data_pos_;
	std::vector<GLVec3> data_col_;

	// list of primitive call (of each kind)
	std::vector<PrimParam> begins_point_;
	std::vector<PrimParam> begins_round_point_;
	std::vector<PrimParam> begins_balls_;
	std::vector<PrimParam> begins_line_;
	std::vector<PrimParam> begins_bold_line_;
	std::vector<PrimParam> begins_face_;
	std::vector<PrimParam>* current_begin_;

	float32 current_size_;
	bool current_aa_;
	bool current_ball_;

public:
	class CGOGN_RENDERING_EXPORT Renderer
	{
		friend class DisplayListDrawer;

		std::unique_ptr<ShaderFlatColorPerVertex::Param> param_cpv_;
		std::unique_ptr<ShaderBoldLineColor::Param> param_bl_;
		std::unique_ptr<ShaderRoundPointColor::Param> param_rp_;
		std::unique_ptr<ShaderPointSpriteColor::Param> param_ps_;
		DisplayListDrawer* drawer_data_;

		Renderer(DisplayListDrawer* dr);

	public:
		~Renderer();

		/**
		 * draw the compiled drawing list
		 * @param projection projection matrix
		 * @param modelview modelview matrix
		 */
		void draw(const GLMat4& projection, const GLMat4& modelview);
	};

	/**
	 * constructor, init all buffers (data and OpenGL) and shader
	 * @Warning need OpenGL context
	 */
	DisplayListDrawer();

	/**
	 * release buffers and shader
	 */
	~DisplayListDrawer();

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(DisplayListDrawer);

	/**
	 * @brief generate a renderer (one per context)
	 * @return pointer on renderer
	 */
	inline std::unique_ptr<Renderer> generate_renderer()
	{
		return std::unique_ptr<Renderer>(new Renderer(this));
	}

	/**
	 * init the data structure
	 */
	void new_list();

	/**
	 * as glBegin, but need a newList call before
	 * @param mode: GL_POINTS, GL_LINES, GL_LINE_LOOP, GL_TRIANGLES,
	 */
	void begin(GLenum mode);

	/**
	 * as glEnd
	 */
	void end();

	/**
	 * finalize the data initialization
	 * drawn is done if newList called with GL_COMPILE_AND_EXECUTE
	 */
	void end_list();

	void vertex3ff(float32 x, float32 y, float32 z);

	/**
	 * use as glVertex3f
	 */
	template <typename SCALAR>
	inline void vertex3f(SCALAR x, SCALAR y, SCALAR z)
	{
		static_assert(std::is_arithmetic<SCALAR>::value, "scalar value only allowed for vertex3");
		vertex3ff(float32(x), float32(y), float32(z));
	}

	void color3ff(float32 r, float32 g, float32 b);

	/**
	 * use as glColor3f
	 */
	template <typename SCALAR>
	inline void color3f(SCALAR x, SCALAR y, SCALAR z)
	{
		static_assert(std::is_arithmetic<SCALAR>::value, "scalar value only allowed for vertex3");
		color3ff(float32(x), float32(y), float32(z));
	}

	inline void vertex3fv(const std::array<float32, 3>& xyz)
	{
		vertex3ff(xyz[0], xyz[1], xyz[2]);
	}

	inline void color3fv(const std::array<float32, 3>& rgb)
	{
		color3ff(rgb[0], rgb[1], rgb[2]);
	}

	template <typename SCALAR>
	inline void vertex3fv(SCALAR* xyz)
	{
		static_assert(std::is_arithmetic<SCALAR>::value, "scalar vector only allowed for vertex3");
		vertex3ff(float32(xyz[0]), float32(xyz[1]), float32(xyz[2]));
	}

	template <typename SCALAR>
	inline void color3fv(SCALAR* rgb)
	{
		static_assert(std::is_arithmetic<SCALAR>::value, "scalar vector only allowed for vertex3");
		color3ff(float32(rgb[0]), float32(rgb[1]), float32(rgb[2]));
	}

	template <typename VEC3>
	inline void vertex3fv(const VEC3& xyz)
	{
		vertex3f(float32(xyz[0]), float32(xyz[1]), float32(xyz[2]));
	}

	template <typename VEC3>
	inline void color3fv(const VEC3& rgb)
	{
		color3f(float32(rgb[0]), float32(rgb[1]), float32(rgb[2]));
	}

	/**
	 * use as glPointSize
	 */
	inline void point_size(float32 ps)
	{
		current_aa_ = false;
		current_size_ = ps;
		current_ball_ = false;
	}

	/**
	 * use as glPointSize with use of anti-aliasing (alpha blending)
	 */
	inline void point_size_aa(float32 ps)
	{
		current_aa_ = true;
		current_size_ = ps;
		current_ball_ = false;
	}

	/**
	 * use as glPointSize for shaded ball rendering
	 */
	inline void ball_size(float32 ps)
	{
		current_ball_ = true;
		current_aa_ = false;
		current_size_ = ps;
	}

	/**
	 * use as glLineWidth
	 */
	inline void line_width(float32 lw)
	{
		current_aa_ = false;
		current_size_ = lw;
	}

	/**
	 * use as glLineWidth with use of anti-aliasing (alpha blending)
	 */
	inline void line_width_aa(float32 lw)
	{
		current_aa_ = true;
		current_size_ = 2.0f * lw;
	}
};

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_DRAWER_H_

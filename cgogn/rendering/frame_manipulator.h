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

#ifndef CGOGN_RENDERING_FRAMEMANIPULATOR_H_
#define CGOGN_RENDERING_FRAMEMANIPULATOR_H_

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/rendering/shaders/frame_manip_drawer.h>
#include <cgogn/rendering/shaders/shader_bold_line.h>
#include <cgogn/rendering/shaders/shader_no_illum.h>
#include <cgogn/rendering/types.h>

namespace cgogn
{

namespace rendering
{

/**
 * @brief The FrameManipulator class
 *
 * Typical usage:
 *
 *  std::unique_ptr<cgogn::rendering::FrameManipulator> frame_manip_;
 *
 * init:
 *   frame_manip_ = std::make_unique<cgogn::rendering::FrameManipulator>();
 *   frame_manip_->set_size(...);
 *   frame_manip_->set_position(...);
 * draw:
 *   frame_manip_->draw(true/false, true/false, proj, view);
 * mousePressEvent:
 *   frame_manip_->pick(event->x(), event->y(),P,Q); // P,Q computed ray
 * mouseReleaseEvent:
 *   frame_manip_->release();
 * mouseMouseEvent:
 *   frame_manip_->drag(true/false, event->x(), event->y());
 *   // for Z plane equation:
 *   Vec3 position;
 *   Vec3 axis_z;
 *   frame_manip_->get_position(position);
 *   frame_manip_->get_axis(cgogn::rendering::FrameManipulator::Zt,axis_z);
 *   float32 d = -(position.dot(axis_z));
 */
class CGOGN_RENDERING_EXPORT FrameManipulator
{

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	std::unique_ptr<VBO> vbo_grid_;
	std::unique_ptr<ShaderNoIllum::Param> param_grid_;

	FrameManipDrawer* fmd_ = FrameManipDrawer::generate();

	static const uint32 nb_grid_ = 10;
	static const uint32 nb_grid_ind_ = 4 * (nb_grid_ + 1);

public:
	enum AXIS : int32
	{
		NONE = 0,
		CENTER,
		Xt,
		Yt,
		Zt,
		Xr,
		Yr,
		Zr,
		Xs,
		Ys,
		Zs,
		Translations,
		Rotations,
		Scales
	};

protected:
	/**
	 * number of segment for circle drawing
	 */
	static const uint32 nb_segments = 64;

	static const float32 ring_half_width;

	/**
	 * locking table
	 */
	bool locked_axis_[11];

	/**
	 * pinking only locking table
	 */
	bool locked_picking_axis_[11];

	/**
	 * VBO for position
	 */
	std::unique_ptr<VBO> vbo_frame_;

	/**
	 * Shader
	 */
	std::unique_ptr<ShaderBoldLine::Param> param_bl_;
	std::unique_ptr<ShaderNoIllum::Param> param_sc_;

	/**
	 * the selectd axis, highlighted
	 */
	uint32 highlighted_;
	bool axis_orientation_;
	GLMat4 rotations_;
	float32 scale_rendering_;
	GLVec3 trans_;
	GLVec3 scale_;
	GLVec3 length_axes_;
	GLVec3 projected_selected_axis_;
	GLVec3 projected_origin_;

	GLMat4 proj_mat_;
	GLMat4 view_mat_;
	GLint viewport_[4];

	// last mouse position
	int beg_X_;
	int beg_Y_;

	inline bool axis_pickable(uint32 a)
	{
		return (!locked_axis_[a]) && (!locked_picking_axis_[a]);
	}

	GLMat4 transfo_render_frame();
	void set_length_axes();
	uint32 pick_frame(const GLVec4& PP, const GLVec4& QQ);
	void store_projection(uint32 ax);
	float32 angle_from_mouse(int x, int y, int dx, int dy);
	float32 distance_from_mouse(int dx, int dy);
	float32 scale_from_mouse(int dx, int dy);
	void translate_in_screen(int dx, int dy);
	void rotate_in_screen(int dx, int dy);

public:
	FrameManipulator();

	/**
	 * set size of frame (for rendering)
	 */
	void set_size(float32 radius);

	/**
	 * @brief set z_plane parameter drawing
	 * @param color
	 * @param xc x position [-1/1] (default 0)
	 * @param yc y position [-1/1] (default 0)
	 * @param r radius (default 1)
	 */
	void z_plane_param(const GLColor& color, float32 xc, float32 yc, float32 r);

	/**
	 * get the size of frame
	 */
	float32 get_size();

	/**
	 * @brief draw the frame and the Z plane
	 * @param frame draw frame or not
	 * @param zplane draw z-plane or not
	 * @param proj projection matrix
	 * @param view model-view matrix
	 */
	void draw(bool frame, bool zplane, const GLMat4& proj, const GLMat4& view);

	/**
	 * @brief try picking the frame
	 * @param x mouse x pos
	 * @param y mouse y pos
	 * @param PP first point of ray picking
	 * @param QQ second point of ray picking
	 */
	template <typename VEC>
	void pick(int x, int y, const VEC& PP, const VEC& QQ);

	/**
	 * @brief release dragging frame
	 */
	inline void release()
	{
		highlighted_ = NONE;
	}

	/**
	 * @brief dragging frame
	 * @param local local frame (true) / screen (false) dragging
	 * @param x mouse position
	 * @param y mouse position
	 */
	void drag(bool local, int x, int y);

	/**
	 * lock an axis (drawing & picking are disabled)
	 * @param axis the axis Xt/Yt/Zt/Xs/Yx/Zs/Xr/Yr/Zr or group Translations/Scales/Rotations)
	 */
	void lock(uint32 axis);

	/**
	 * unlock an axis
	 * @param axis the axis Xt/Yt/Zt/Xs/Yx/Zs/Xr/Yr/Zr or group Translations/Scales/Rotations)
	 */
	void unlock(uint32 axis);

	/**
	 * is an axis locked
	 * @param axis the axis to test
	 */
	bool locked(uint32 axis);

	/**
	 * lock an axis only for pinking
	 */
	void lock_picking(uint32 axis);

	/**
	 * unlock an axis (only for pinking)
	 */
	void unlock_picking(uint32 axis);

	/**
	 * is an axis locked (only for pinking)
	 */
	bool locked_picking(uint32 axis);

	/**
	 * higlight an axis (change width rendering).
	 * To unhighlight, just highlight NONE or highlight a already highlighted  axis
	 */
	void highlight(uint32 axis);

	/**
	 * rotate the frame around one of its axis
	 */
	void rotate(uint32 axis, float32 angle);

	/**
	 * translate the frame around one of its axis
	 * @param axis
	 * @param x ratio of frame radius
	 */
	void translate(uint32 axis, float32 x);

	/**
	 * scale the frame in direction of one axis
	 * @param axis (Xs/Ys/Zs/CENTER)
	 * @param sc scale factor to apply on
	 */
	void scale(uint32 axis, float32 sc);

	/**
	 * get the matrix transformation
	 */
	GLMat4 transfo();

	/**
	 * set the position of frame
	 * @param P the position of origin
	 */
	template <typename VEC3>
	void set_position(const VEC3& P);

	inline GLVec3 get_position()
	{
		return trans_;
	}

	// get position in a non-QVector3 vector
	template <typename VEC3>
	void get_position(VEC3& pos);

	/**
	 * @brief get an axis
	 * @param ax (Xr,Yr,Zr)
	 * @return the axis
	 */
	GLVec3 get_axis(uint32 ax);

	// get axis in a non-QVector3 vector
	template <typename VEC3>
	void get_axis(uint32 ax, VEC3& pos);

	/**
	 * set the scale of frame
	 * @param P the vector of scale factors
	 */
	void set_scale(const GLVec3& S);

	/**
	 * set the orientation of frame (Z is deduced)
	 * @param X the vector X of frame
	 * @param Y the vector Y of frame
	 * @return return false if parameters are not unit orthogonal vectors
	 */
	bool set_orientation(const GLVec3& X, const GLVec3& Y);

	/**
	 * set transformation matrix
	 */
	void set_transformation(const GLMat4& transfo);

	inline static bool rotation_axis(uint32 axis)
	{
		return (axis >= Xr) && (axis <= Zr);
	}

	inline static bool translation_axis(uint32 axis)
	{
		return (axis >= Xt) && (axis <= Zt);
	}

	inline static bool scale_axis(uint32 axis)
	{
		return ((axis >= Xs) && (axis <= Zs)) || (axis == CENTER);
	}
};

// template implementation

template <typename VEC>
void FrameManipulator::pick(int x, int y, const VEC& PP, const VEC& QQ)
{
	beg_X_ = x;
	beg_Y_ = y;

	GLVec4 P(PP[0], PP[1], PP[2], 1.0);
	GLVec4 Q(QQ[0], QQ[1], QQ[2], 1.0);
	highlighted_ = pick_frame(P, Q);

	if (highlighted_ != NONE)
		store_projection(highlighted_);
}

template <typename VEC3>
void FrameManipulator::set_position(const VEC3& pos)
{
	trans_[0] = pos[0];
	trans_[1] = pos[1];
	trans_[2] = pos[2];
}

template <typename VEC3>
void FrameManipulator::get_position(VEC3& pos)
{
	pos[0] = trans_[0];
	pos[1] = trans_[1];
	pos[2] = trans_[2];
}

template <typename VEC3>
void FrameManipulator::get_axis(uint32 ax, VEC3& axis)
{
	GLVec3 A = get_axis(ax);
	axis[0] = A[0];
	axis[1] = A[1];
	axis[2] = A[2];
}

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_FRAMEMANIPULATOR_H_

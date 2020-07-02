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

#include <cgogn/rendering/frame_manipulator.h>
#include <cgogn/rendering/vbo_update.h>

#include <cgogn/geometry/functions/distance.h>
#include <cgogn/geometry/functions/intersection.h>

// #include <cgogn/geometry/types/vector_traits.h>

#include <cmath>

namespace cgogn
{

namespace rendering
{

const float32 FrameManipulator::ring_half_width = 0.08f;

FrameManipulator::FrameManipulator()
	: highlighted_(NONE), scale_rendering_(1.0f), trans_(0.0f, 0.0f, 0.0f), scale_(1.0f, 1.0f, 1.0f)
{
	fmd_ = FrameManipDrawer::generate();
	rotations_.setIdentity();
	for (uint32 i = 0; i < 11; ++i)
	{
		locked_axis_[i] = false;
		locked_picking_axis_[i] = false;
	}
	set_length_axes();
}

void FrameManipulator::set_size(float32 radius)
{
	if (scale_rendering_ >= 0.0f)
		scale_rendering_ = radius;
}

float32 FrameManipulator::get_size()
{
	return scale_rendering_;
}

void FrameManipulator::draw(bool frame, bool zplane, const GLMat4& proj, const GLMat4& view)
{
	fmd_->set_axis_selected(highlighted_ - Xt);
	fmd_->set_ring_selected(highlighted_ - Xr);
	GLMat4 fr = transfo_render_frame();

	if (frame)
	{
		fmd_->draw_transla(proj, view, fr);
		fmd_->draw_rota(proj, view, fr);
	}
	if (zplane)
		fmd_->draw_grid(proj, view, fr);

	proj_mat_ = proj;
	view_mat_ = view;
	glGetIntegerv(GL_VIEWPORT, viewport_);
}

void FrameManipulator::highlight(uint32 axis)
{
	if (highlighted_ == axis)
		highlighted_ = NONE;
	else
		highlighted_ = axis;
}

uint32 FrameManipulator::pick_frame(const GLVec4& PP, const GLVec4& QQ)
{
	// ray inverse transfo
	GLMat4 invtr = transfo_render_frame().inverse();
	GLVec4 tP = invtr * PP;
	GLVec4 tQ = invtr * QQ;

	Vec3 P(tP[0] / tP[3], tP[1] / tP[3], tP[2] / tP[3]);
	Vec3 Q(tQ[0] / tQ[3], tQ[1] / tQ[3], tQ[2] / tQ[3]);
	Vec3 V = Q - P;

	// origin of frame
	Vec3 origin(0.0, 0.0, 0.0);

	// intersection possible between line and frame (10% margin)?
	float32 dist2 = float32(cgogn::geometry::squared_distance_line_point(P, Q, origin));

	float32 distMax = std::max(length_axes_[0], std::max(length_axes_[1], length_axes_[2]));
	distMax *= 3.6f;
	distMax = std::max(distMax, 1.0f + ring_half_width);

	if (dist2 > distMax * distMax)
		return NONE;

	// click on center
	if (dist2 < 0.02f * 0.02f)
	{
		if (axis_pickable(CENTER))
			return CENTER;
		else
			return NONE;
	}

	Scalar dist_target[9];
	Scalar dist_cam[9];

	for (uint32 i = 0; i < 9; ++i)
	{
		dist_target[i] = 2.0f;
		dist_cam[i] = std::numeric_limits<Scalar>::max();
	}

	// Circles:
	// plane X=0
	Vec3 Qx;
	bool inter = cgogn::geometry::intersection_line_plane(P, V, origin, Vec3(1.0, 0.0, 0.0), &Qx);

	if (axis_pickable(Xr))
	{
		if (inter)
			dist_target[3] = Qx.norm() - 1.0;

		if (std::abs(dist_target[3]) < ring_half_width)
			dist_cam[3] = (P - Qx).squaredNorm();
	}

	// plane Y=0
	Vec3 Qy;
	inter = cgogn::geometry::intersection_line_plane(P, V, origin, Vec3(0.0, 1.0, 0.0), &Qy);

	if (axis_pickable(Yr))
	{
		if (inter)
			dist_target[4] = Qy.norm() - 1.0;

		if (std::abs(dist_target[4]) < ring_half_width)
			dist_cam[4] = (P - Qy).squaredNorm();
	}

	// plane Z=0
	Vec3 Qz;
	inter = cgogn::geometry::intersection_line_plane(P, V, origin, Vec3(0.0, 0.0, 1.0), &Qz);

	if (axis_pickable(Zr))
	{
		if (inter)
			dist_target[5] = Qz.norm() - 1.0;

		if (std::abs(dist_target[5]) < ring_half_width)
			dist_cam[5] = (P - Qz).squaredNorm();
	}

	// axes:

	if (axis_pickable(Xt) || axis_pickable(Xs))
	{
		Vec3 PX(3.6 * length_axes_[0], 0.0, 0.0);
		dist_target[0] = std::sqrt(cgogn::geometry::squared_distance_line_seg(P, V, V.squaredNorm(), origin, PX));
		if (std::abs(dist_target[0]) < 0.02)
		{
			if (axis_pickable(Xt) && !axis_pickable(Xs))
				dist_cam[0] = (P - PX).squaredNorm();
			else
			{
				if (Qz.norm() > length_axes_[0])
					dist_cam[0] = (P - PX).squaredNorm();
				else
					dist_cam[6] = P.squaredNorm();
			}
		}
	}

	if (axis_pickable(Yt) || axis_pickable(Ys))
	{
		Vec3 PY(0.0, 3.6 * length_axes_[1], 0.0);
		dist_target[1] = std::sqrt(cgogn::geometry::squared_distance_line_seg(P, V, V.squaredNorm(), origin, PY));
		if (std::abs(dist_target[1]) < 0.02)
		{
			if (axis_pickable(Yt) && !axis_pickable(Ys))
				dist_cam[1] = (P - PY).squaredNorm();
			else
			{
				if (Qz.norm() > length_axes_[1])
					dist_cam[1] = (P - PY).squaredNorm();
				else
					dist_cam[7] = P.squaredNorm();
			}
		}
	}

	if (axis_pickable(Zt) || axis_pickable(Zs))
	{
		Vec3 PZ(0.0, 0.0, 3.6 * length_axes_[2]);
		dist_target[2] = std::sqrt(cgogn::geometry::squared_distance_line_seg(P, V, V.squaredNorm(), origin, PZ));
		if (std::abs(dist_target[2]) < 0.02)
		{
			if (axis_pickable(Zt) && !axis_pickable(Zs))
				dist_cam[2] = (P - PZ).squaredNorm();
			else
			{
				if (Qx.norm() > length_axes_[2])
					dist_cam[2] = (P - PZ).squaredNorm();
				else
					dist_cam[8] = P.squaredNorm();
			}
		}
	}

	// find min dist_cam value;
	uint32 min_index = 0;
	Scalar min_val = dist_cam[0];
	for (uint32 i = 1; i < 9; ++i)
	{
		if (dist_cam[i] < min_val)
		{
			min_val = dist_cam[i];
			min_index = i;
		}
	}

	if (min_val < std::numeric_limits<Scalar>::max())
	{
		if (axis_pickable(Xt + min_index))
			return Xt + min_index;
	}

	return NONE;
}

void FrameManipulator::rotate(uint32 axis, float32 angle)
{
	// create axis
	GLVec3 ax(0, 0, 0);
	ax[axis - Xr] = 1.0f;

	Eigen::Transform<float, 3, Eigen::Affine> tr(Eigen::AngleAxisf(angle, ax));
	rotations_ = rotations_ * tr.matrix();
}

void FrameManipulator::translate(uint32 axis, float32 x)
{
	trans_ += x * scale_rendering_ * rotations_.block<3, 1>(0, axis - Xt);
}

void FrameManipulator::set_length_axes()
{
	float32 avgScale = 0.75f / (scale_[0] + scale_[1] + scale_[2]);

	length_axes_ = avgScale * scale_;
}

void FrameManipulator::scale(uint32 axis, float32 sc)
{
	if (axis == CENTER)
	{
		scale_[0] *= sc;
		scale_[1] *= sc;
		scale_[2] *= sc;
	}
	else
		scale_[axis - Xs] *= sc;

	set_length_axes();
}

GLMat4 FrameManipulator::transfo_render_frame()
{
	GLMat4 tr = GLMat4::Identity();
	tr.block<3, 1>(0, 3) = trans_;

	tr *= rotations_;

	float32 avgScale = (scale_[0] + scale_[1] + scale_[2]) / 3.0f * scale_rendering_;
	Eigen::Transform<float, 3, Eigen::Affine> sc(Eigen::Scaling(avgScale));
	tr *= sc.matrix();
	return tr;
}

GLMat4 FrameManipulator::transfo()
{
	GLMat4 tr = GLMat4::Identity();
	tr.block<3, 1>(0, 3) = trans_;
	tr *= rotations_;
	Eigen::Transform<float, 3, Eigen::Affine> sc(Eigen::Scaling(scale_));
	tr *= sc.matrix();
	return tr;
}

void FrameManipulator::set_scale(const GLVec3& S)
{
	scale_ = S;
	set_length_axes();
}

bool FrameManipulator::set_orientation(const GLVec3& X, const GLVec3& Y)
{
	GLVec3 Z = X.cross(Y);

	if ((X.norm() != 1.0f) || (Y.norm() != 1.0f) || (Z.norm() != 1.0f))
		return false;

	rotations_.block<3, 1>(0, 0) = X;
	rotations_.block<3, 1>(0, 1) = Y;
	rotations_.block<3, 1>(0, 2) = Z;

	return true;
}

void FrameManipulator::set_transformation(const GLMat4&)
{
	// TODO E.S.: parameter is not used. It seems wrong.
	set_position(rotations_.block<3, 1>(0, 3).eval());

	//	col = rotations_.column(0);
	//	QVector3D Rx(	col[0], col[1], col[2]);
	//	col = rotations_.column(1);
	//	QVector3D Ry(	col[0], col[1], col[2]);
	//	col = rotations_.column(2);
	//	QVector3D Rz(	col[0], col[1], col[2]);

	//	set_scale(QVector3D(float32(Rx.length()), float32(Ry.length()), float32(Rz.length())));

	//	col[3] = 0.0f;
	//	col[0] = Rx[0];
	//	col[1] = Rx[1];
	//	col[2] = Rx[2];
	//	rotations_.setColumn(0,col);

	//	col[0] = Ry[0];
	//	col[1] = Ry[1];
	//	col[2] = Ry[2];
	//	rotations_.setColumn(0,col);

	//	col[0] = Rz[0];
	//	col[1] = Rz[1];
	//	col[2] = Rz[2];
	//	rotations_.setColumn(0,col);
}

void FrameManipulator::lock(uint32 axis)
{
	assert(axis <= Scales);
	switch (axis)
	{
	case Translations:
		locked_axis_[Xt] = true;
		locked_axis_[Yt] = true;
		locked_axis_[Zt] = true;
		break;
	case Rotations:
		locked_axis_[Xr] = true;
		locked_axis_[Yr] = true;
		locked_axis_[Zr] = true;
		break;
	case Scales:
		locked_axis_[Xs] = true;
		locked_axis_[Ys] = true;
		locked_axis_[Zs] = true;
		break;
	default:
		locked_axis_[axis] = true;
		break;
	}
	set_length_axes();
}

void FrameManipulator::unlock(uint32 axis)
{
	assert(axis <= Scales);
	switch (axis)
	{
	case Translations:
		locked_axis_[Xt] = false;
		locked_axis_[Yt] = false;
		locked_axis_[Zt] = false;
		break;
	case Rotations:
		locked_axis_[Xr] = false;
		locked_axis_[Yr] = false;
		locked_axis_[Zr] = false;
		break;
	case Scales:
		locked_axis_[Xs] = false;
		locked_axis_[Ys] = false;
		locked_axis_[Zs] = false;
		break;
	default:
		locked_axis_[axis] = false;
		break;
	}
	set_length_axes();
}

bool FrameManipulator::locked(uint32 axis)
{
	assert(axis <= Zs);
	return locked_axis_[axis];
}

void FrameManipulator::lock_picking(uint32 axis)
{
	assert(axis <= Scales);
	switch (axis)
	{
	case Translations:
		locked_picking_axis_[Xt] = true;
		locked_picking_axis_[Yt] = true;
		locked_picking_axis_[Zt] = true;
		break;
	case Rotations:
		locked_picking_axis_[Xr] = true;
		locked_picking_axis_[Yr] = true;
		locked_picking_axis_[Zr] = true;
		break;
	case Scales:
		locked_picking_axis_[Xs] = true;
		locked_picking_axis_[Ys] = true;
		locked_picking_axis_[Zs] = true;
		break;
	default:
		locked_picking_axis_[axis] = true;
		break;
	}
	set_length_axes();
}

void FrameManipulator::unlock_picking(uint32 axis)
{
	assert(axis <= Scales);
	switch (axis)
	{
	case Translations:
		locked_picking_axis_[Xt] = false;
		locked_picking_axis_[Yt] = false;
		locked_picking_axis_[Zt] = false;
		break;
	case Rotations:
		locked_picking_axis_[Xr] = false;
		locked_picking_axis_[Yr] = false;
		locked_picking_axis_[Zr] = false;
		break;
	case Scales:
		locked_picking_axis_[Xs] = false;
		locked_picking_axis_[Ys] = false;
		locked_picking_axis_[Zs] = false;
		break;
	default:
		locked_picking_axis_[axis] = false;
		break;
	}
	set_length_axes();
}

bool FrameManipulator::locked_picking(uint32 axis)
{
	return locked_picking_axis_[axis];
}

GLVec3 FrameManipulator::get_axis(uint32 ax)
{
	uint32 i = (ax - Xt) % 3;
	return rotations_.block<3, 1>(0, i);
}

void FrameManipulator::store_projection(uint32 ax)
{
	GLMat4 mat = proj_mat_ * view_mat_;
	GLVec4 tr(trans_.x(), trans_.y(), trans_.z(), 1);
	GLVec4 Or4 = (mat * tr);
	GLVec3 Or = Or4.head<3>() / Or4[3];
	projected_origin_ = Or * 0.5f + GLVec3(0.5f, 0.5f, 0.5f);
	projected_origin_[0] = projected_origin_[0] * float32(viewport_[2]) + float32(viewport_[0]);
	projected_origin_[1] = projected_origin_[1] * float32(viewport_[3]) + float32(viewport_[1]);

	if (ax > CENTER)
	{
		GLVec3 A = get_axis(ax);
		if ((ax >= Xr) && (ax <= Zr))
		{
			// compute screen orientation
			GLVec3 V = (view_mat_.block<3, 3>(0, 0) * A);
			axis_orientation_ = V[2] > 0;
		}

		A = A + trans_;
		GLVec4 A4(A.x(), A.y(), A.z(), 1);
		A4 = (mat * A4);
		projected_selected_axis_ = A4.head<3>() / A4[3];
		projected_selected_axis_ = projected_selected_axis_ * 0.5f + GLVec3(0.5f, 0.5f, 0.5f);
		projected_selected_axis_[0] = projected_selected_axis_[0] * float32(viewport_[2]) + float32(viewport_[0]);
		projected_selected_axis_[1] = projected_selected_axis_[1] * float32(viewport_[3]) + float32(viewport_[1]);
		projected_selected_axis_ -= projected_origin_;
	}
}

float32 FrameManipulator::angle_from_mouse(int x, int y, int dx, int dy)
{
	Vec3 Vo(float32(x) - projected_origin_[0], float32(viewport_[3] - y) - projected_origin_[1], 0.0f);

	Vec3 dV(float32(dx), float32(dy), 0.0f);

	Vo.normalize();
	dV.normalize();
	Vec3 W = Vo.cross(dV);
	return (axis_orientation_ ? W[2] : -W[2]) / 100.0f;
}

float32 FrameManipulator::distance_from_mouse(int dx, int dy)
{
	Vec3 dV(float32(dx), float32(dy), 0.0f);
	Vec3 ax(projected_selected_axis_[0], projected_selected_axis_[1], projected_selected_axis_[2]);
	float32 tr = float32(dV.dot(ax));

	if (tr > 0)
		tr = float32(dV.norm() / 100.0f);
	else
		tr = float32(dV.norm() / -100.0f);

	return tr;
}

float32 FrameManipulator::scale_from_mouse(int dx, int dy)
{
	if (abs(dx) > abs(dy))
	{
		if (dx > 0)
			return 1.01f;
		return 0.99f;
	}
	else
	{
		if (dy > 0)
			return 1.01f;
		return 0.99f;
	}
}

void FrameManipulator::translate_in_screen(int dx, int dy)
{
	// unproject origin shifted
	GLMat4 inv_mat = (proj_mat_ * view_mat_).inverse();
	GLVec4 P(projected_origin_[0] + float32(dx), projected_origin_[1] + float32(dy), projected_origin_[2], 1.0f);
	P[0] = (P[0] - float32(viewport_[0])) / float32(viewport_[2]);
	P[1] = (P[1] - float32(viewport_[1])) / float32(viewport_[3]);
	P = P * 2.0f - GLVec4(1.0f, 1.0f, 1.0f, 1.0f);

	P = inv_mat * P;
	P /= P[3];

	// and store it
	trans_[0] = P[0];
	trans_[1] = P[1];
	trans_[2] = P[2];
	store_projection(NONE);
}

void FrameManipulator::rotate_in_screen(int dx, int dy)
{
	GLMat4 inv_mat = (proj_mat_ * view_mat_).inverse();
	GLVec4 P(dx, -dy, 0.0f, 1.0f);
	P = inv_mat * P;
	P /= P[3];

	GLVec3 axis_rot(P[1], P[0], P[2]);
	axis_rot.normalize();

	float32 angle = std::sqrt(float32(dx * dx + dy * dy)) / 200.0f;
	Eigen::Transform<float, 3, Eigen::Affine> tr(Eigen::AngleAxisf(angle, axis_rot));
	rotations_ = tr * rotations_;
}

void FrameManipulator::drag(bool local, int x, int y)
{
	// rotation selected ?
	if (rotation_axis(highlighted_))
	{
		if (local)
		{
			float angle = angle_from_mouse(x, y, x - beg_X_, beg_Y_ - y);
			rotate(highlighted_, angle);
		}
		else
			rotate_in_screen(x - beg_X_, beg_Y_ - y);
	}
	// translation selected
	else if (translation_axis(highlighted_))
	{
		if (local)
		{
			float dist = distance_from_mouse(x - beg_X_, beg_Y_ - y);
			translate(highlighted_, dist);
		}
		else
			translate_in_screen(x - beg_X_, beg_Y_ - y);
	}
	// scale selected
	else if (scale_axis(highlighted_))
	{
		float sc = scale_from_mouse(x - beg_X_, beg_Y_ - y);
		scale(highlighted_, sc);
	}

	beg_X_ = x;
	beg_Y_ = y;
}

} // namespace rendering

} // namespace cgogn

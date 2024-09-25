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

#ifndef CGOGN_RENDERING_TYPES_H_
#define CGOGN_RENDERING_TYPES_H_

#include <cgogn/core/utils/numerics.h>
#include <cgogn/rendering/cgogn_rendering_export.h>

#include <Eigen/Dense>

#include <GL/gl3w.h>
#include <iostream>
#include <map>
#include <string>

namespace cgogn
{

namespace rendering
{

using GLVec2 = Eigen::Vector2f;
using GLVec3 = Eigen::Vector3f;
using GLVec4 = Eigen::Vector4f;

template <typename T1, typename T2>
inline GLVec2 construct_GLVec4(T1 x, T2 y)
{
	return GLVec2(float32(x), float32(y));
}

template <typename T1, typename T2, typename T3>
inline GLVec4 construct_GLVec4(T1 x, T2 y, T3 z)
{
	return GLVec4(float32(x), float32(y), float32(z));
}

template <typename T1, typename T2, typename T3, typename T4>
inline GLVec4 construct_GLVec4(T1 x, T2 y, T3 z, T4 w)
{
	return GLVec4(float32(x), float32(y), float32(z), float32(w));
}

using GLVec2d = Eigen::Vector2d;
using GLVec3d = Eigen::Vector3d;
using GLVec4d = Eigen::Vector4d;

using GLMat3 = Eigen::Matrix3f;
using GLMat4 = Eigen::Matrix4f;

using GLMat3d = Eigen::Matrix3d;
using GLMat4d = Eigen::Matrix4d;

using GLColor = Eigen::Vector4f;

using Transfo3d = Eigen::Affine3d;

using Transfo3f = Eigen::Affine3f;

// inline GLColor col4i(uint8 R, uint8 G, uint8 B, uint8 A)
// {
// 	return GLColor(R/255.0f, G/255.0f, B/255.0f, A/255.0f);
// }

// inline GLColor col3i(uint8 R, uint8 G, uint8 B)
// {
// 	return GLColor(R/255.0f, G/255.0f, B/255.0f, 1.0f);
// }

// inline GLColor col1i(uint8 R)
// {
// 	return GLColor(R/255.0f, R/255.0f, R/255.0f, 1.0f);
// }

static std::map<GLenum, std::string> GL_ERRORS_NAMES = {
	{GL_INVALID_ENUM, "GL_INVALID_ENUM"},
	{GL_INVALID_VALUE, "GL_INVALID_VALUE"},
	{GL_INVALID_OPERATION, "GL_INVALID_OPERATION"},
	{GL_INVALID_FRAMEBUFFER_OPERATION, "GL_INVALID_FRAMEBUFFER_OPERATION"},
	{GL_OUT_OF_MEMORY, "GL_OUT_OF_MEMORY"},
	{GL_STACK_UNDERFLOW, "GL_STACK_UNDERFLOW"},
	{GL_STACK_OVERFLOW, "GL_STACK_OVERFLOW"}};

#ifdef CGOGN_GL43_DEBUG_MODE
inline void gl_debug_name(GLenum type, GLuint id, const std::string& name)
{
	glObjectLabel(type, id, GLsizei(name.length()), name.c_str());
}
#else
inline void gl_debug_name(GLenum, GLuint, const std::string&)
{
}
#endif

} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_TYPES_H_

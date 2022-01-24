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

#ifndef CGOGN_RENDERING_SHADERS_FUNCTION_COLORMAPS_H_
#define CGOGN_RENDERING_SHADERS_FUNCTION_COLORMAPS_H_

#include <string>

namespace cgogn
{

namespace rendering
{

enum ColorMap
{
	BWR = 0,
	CWR,
	BCGYR,
	BGR
};

namespace shader_function
{

struct ColorMap
{
	static const std::string source;

	struct Uniforms
	{
		int color_map_;
		int expansion_;
		float min_value_;
		float max_value_;

		inline Uniforms() : color_map_(BWR), expansion_(0), min_value_(0), max_value_(1)
		{
		}
	};

	static const char* uniform_names[4];
};

} // namespace shader_function

} // namespace rendering

} // namespace cgogn

#endif

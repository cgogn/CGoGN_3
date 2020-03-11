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

namespace shader_funcion
{

inline std::string color_maps_shader_source()
{
	return std::string(R"(
	const float M_PI = 3.1415926535897932384626433832795;

	uniform int color_map;
	uniform int expansion;
	uniform float min_value;
	uniform float max_value;

	float scale_and_clamp_to_0_1(float x, float min, float max)
	{
		// smooth_step ?? (hermite instead of linear)
		// smooth_step(x,min,max);
		float v = (x - min) / (max - min);
		return clamp(v, 0.0, 1.0);
	}

	float scale_expand_within_0_1(float x, int n)
	{
		for (int i = 1; i <= n; i++)
			x = (1.0 - cos(M_PI * x)) / 2.0;
		for (int i = -1; i >= n; i--)
			x = acos(1.0 - 2.0 * x) / M_PI;
		return x;
	}

	vec3 color_map_blue_white_red(float x)
	{
		float x2 = 2.0 * x;
		switch(int(floor(max(0.0,x2+1.0))))
		{
			case 0: return vec3(0.0, 0.0, 1.0) ;
			case 1: return vec3(x2, x2 , 1.0);
			case 2: return vec3(1.0, 2.0 - x2, 2.0 - x2);
		}
		return vec3(1.0, 0.0, 0.0) ;
	}

	vec3 color_map_cyan_white_red(float x)
	{
		float x2 = 2.0 * x;
		switch(int(floor(max(0.0,x2+1.0))))
		{
			case 0: return vec3(0.0, 0.0, 1.0) ;
			case 1: return vec3(x2, 1.0 , 1.0);
			case 2: return vec3(1.0, 2.0 - x2, 2.0 - x2);
		}
		return vec3(1.0, 0.0, 0.0) ;
	}

	vec3 color_map_BCGYR(float x)
	{
		float x4 = 4.0 * x;
		switch(int(floor(max(0.0,x4+1.0))))
		{
			case 0: return vec3(0, 0, 1) ;
			case 1: return vec3(0.0, x4, 1.0);
			case 2: return vec3(0.0, 1.0 , 2.0 - x4);
			case 3: return vec3(x4 - 2.0, 1.0, 0.0);
			case 4: return vec3(1.0, 4.0 - x4, 0.0);
		}
		return vec3(1, 0, 0) ;
	}

	vec3 color_map_blue_green_red(float x)
	{
		float x2 = 2.0 * x;
		switch(int(floor(max(0.0,x2+1.0))))
		{
			case 0: return vec3(0.0, 0.0, 1.0) ;
			case 1: return vec3(0.0, 2.0 * x, 1.0 - 2.0 * x);
			case 2: return vec3(2.0 * x - 1.0, 2.0 - 2.0 * x, 0.0);
		}
		return vec3(1.0, 0.0, 0.0) ;
	}

	vec3 scalar2color(float scalar)
	{
		float value = scale_expand_within_0_1(scale_and_clamp_to_0_1(
					scalar,
					min_value,
					max_value),
				expansion);
		switch(color_map)
		{
			case 0 : return color_map_blue_white_red(value);
			case 1 : return color_map_cyan_white_red(value);
			case 2 : return color_map_BCGYR(value);
			case 3 : return color_map_blue_green_red(value);
		}
		return vec3(1,1,1);
	}

	)");
}
} // namespace shader_funcion

} // namespace rendering
} // namespace cgogn

#endif

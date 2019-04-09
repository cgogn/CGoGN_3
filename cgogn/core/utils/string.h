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

#ifndef CGOGN_CORE_UTILS_STRING_H_
#define CGOGN_CORE_UTILS_STRING_H_

#include <cgogn/core/cgogn_core_export.h>

#include <string>
#include <locale>
#include <iostream>

namespace cgogn
{

CGOGN_CORE_EXPORT std::string trim_left(const std::string& str);

CGOGN_CORE_EXPORT std::string trim_right(const std::string& str);

CGOGN_CORE_EXPORT std::string trim(const std::string& str);

CGOGN_CORE_EXPORT std::string to_upper(const std::string& str);

CGOGN_CORE_EXPORT std::string to_lower(const std::string& str);

CGOGN_CORE_EXPORT std::string extension(const std::string& str);

CGOGN_CORE_EXPORT std::string remove_extension(const std::string& str);

CGOGN_CORE_EXPORT bool i_equals(const std::string& str1, const std::string& str2);

} // namespace cgogn

#endif // CGOGN_CORE_UTILS_STRING_H_

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


#include <cgogn/core/utils/string.h>
#include <algorithm>

namespace cgogn
{

CGOGN_CORE_EXPORT std::string trim_left(const std::string& str)
{
	const std::string pattern = " \f\n\r\t\v";
	return str.substr(str.find_first_not_of(pattern));
}

//
//Right trim
//
CGOGN_CORE_EXPORT std::string trim_right(const std::string& str)
{
	const std::string pattern = " \f\n\r\t\v";
	return str.substr(0,str.find_last_not_of(pattern) + 1);
}

//
//Left and Right trim
//
CGOGN_CORE_EXPORT std::string trim(const std::string& str)
{
	return trim_left(trim_right(str));
}

CGOGN_CORE_EXPORT std::string to_upper(const std::string& str)
{
	const std::locale locale;
	std::string res(str);
	for (auto& c : res)
		c = std::string::value_type(std::toupper(c,locale));
	return res;
}

CGOGN_CORE_EXPORT std::string to_lower(const std::string& str)
{
	const std::locale locale;
	std::string res(str);
	for (auto& c : res)
		c = std::string::value_type(std::tolower(c,locale));
	return res;
}

std::string filename_from_path(const std::string& s)
{
	using char_t = std::string::value_type;
	char_t sep = '/';

#ifdef _WIN32
	sep = '\\';
#endif

	const std::size_t i = s.rfind(sep, s.length());
	if (i != std::string::npos)
		return(s.substr(i+1, s.length() - i));

	return("");
}

CGOGN_CORE_EXPORT std::string extension(const std::string& str)
{
	const std::size_t dot = str.rfind('.');
	if (dot == std::string::npos || dot == str.size() - 1u)
		return std::string();
	return str.substr(dot + 1u);
}

CGOGN_CORE_EXPORT std::string remove_extension(const std::string& str)
{
	const std::size_t dot = str.rfind('.');
	if (dot == std::string::npos || dot == str.size() - 1u)
		return str;
	else
		return str.substr(0, dot);
}

CGOGN_CORE_EXPORT bool i_equals(const std::string& str1, const std::string& str2)
{
	using char_t = std::string::value_type;

	if (str1.size() != str2.size())
		return false;

	const std::locale locale;
	return std::equal(str1.begin(), str1.end(), str2.begin(), [&locale] (char_t c1, char_t c2)
	{
		return std::tolower(c1, locale) == std::tolower(c2, locale);
	});
}

} // namespace cgogn

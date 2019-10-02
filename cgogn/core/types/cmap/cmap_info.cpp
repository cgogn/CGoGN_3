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

#include <cgogn/core/types/cmap/cmap_info.h>

#include <iostream>
#include <iomanip>

namespace cgogn
{

void dump_map(const CMapBase& m)
{
    m.foreach_dart([&] (Dart d) -> bool
    {
        std::cout << "index: " << std::setw(5) << d.index << " / ";
        for (auto& r : m.relations_)
        {
            std::cout << r->name() << ": " << (*r)[d.index] << " / ";
        }
        std::cout << " boundary: " << std::boolalpha << m.is_boundary(d) << std::endl;
        return true;
    });
}

}

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

#ifndef CGOGN_CORE_GMAP_GMAP_BASE_H_
#define CGOGN_CORE_GMAP_GMAP_BASE_H_

#include <cgogn/core/types/cmap/map_base.h>

namespace cgogn
{

struct CGOGN_CORE_EXPORT GMapBase : public MapBase
{
};

} // namespace cgogn

#include <cgogn/core/types/cmap/gmap/gmap_orbit_traversals.hpp>
#include <cgogn/core/types/cmap/cmap_base_index.hpp>
#include <cgogn/core/types/cmap/cmap_base_attribute.hpp>
#include <cgogn/core/types/cmap/gmap/gmap_local_traversals.hpp>
#include <cgogn/core/types/cmap/cmap_base_global_traversals.hpp>
#include <cgogn/core/types/cmap/cmap_base_functions.hpp>

#endif // CGOGN_CORE_CMAP_CMAP_BASE_H_

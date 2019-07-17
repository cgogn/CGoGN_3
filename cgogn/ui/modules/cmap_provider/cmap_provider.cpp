/*******************************************************************************
* CGoGN                                                                        *
* Copyright (C) 2019, IGG Group, ICube, University of Strasbourg, France       *
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

#include <cgogn/ui/modules/cmap_provider/cmap_provider.h>

#include <cgogn/io/surface_import.h>
#include <cgogn/core/utils/string.h>

#include <imgui/imgui.h>

namespace cgogn
{

namespace ui
{

CMapProvider::CMapProvider(const App& app) : Module(app, "CMapProvider")
{}

CMapProvider::~CMapProvider()
{}

CMap2* CMapProvider::import_surface_from_file(const std::string& filename)
{
	const auto [it, inserted] = cmap2_.emplace(filename_from_path(filename), std::make_unique<CMap2>());
	CMap2* mesh = (*it).second.get();
	cgogn::io::import_OFF(*mesh, filename);
	return mesh;
}

void CMapProvider::interface()
{}

} // namespace ui

} // namespace cgogn

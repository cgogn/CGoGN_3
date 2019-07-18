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

#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/io/surface_import.h>
#include <cgogn/core/utils/string.h>

#include <imgui/imgui.h>

namespace cgogn
{

namespace ui
{

MeshProvider::MeshProvider(const App& app) : Module(app, "MeshProvider")
{}

MeshProvider::~MeshProvider()
{}

MeshProvider::Mesh* MeshProvider::import_surface_from_file(const std::string& filename)
{
	const auto [it, inserted] = meshes_.emplace(filename_from_path(filename), std::make_unique<Mesh>());
	Mesh* m = it->second.get();
	cgogn::io::import_OFF(*m, filename);
	mesh_data_.emplace(m, MeshData(m));
	return m;
}

MeshData* MeshProvider::mesh_data(const Mesh* m)
{
	return &mesh_data_[m];
	// if (auto it = mesh_data_.find(m); it != mesh_data_.end())
	// 	return &(it->second);
	// else
	// 	return nullptr;
}

void MeshProvider::interface()
{}

} // namespace ui

} // namespace cgogn

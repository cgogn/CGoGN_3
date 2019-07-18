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

#ifndef CGOGN_MODULE_MESH_PROVIDER_MESH_DATA_H_
#define CGOGN_MODULE_MESH_PROVIDER_MESH_DATA_H_

#include <cgogn/ui/modules/mesh_provider/cgogn_module_mesh_provider_export.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/rendering/mesh_render.h>
#include <cgogn/rendering/vbo_update.h>

namespace cgogn
{

namespace ui
{

struct MeshData
{
    using Mesh = CMap2;

    template <typename T>
    using Attribute = typename mesh_traits<Mesh>::Attribute<T>;

    using Vertex = typename mesh_traits<Mesh>::Vertex;

    using Vec3 = geometry::Vec3;
    using Vec2 = geometry::Vec2;
    using Scalar = geometry::Scalar;

	MeshData() : mesh_(nullptr)
	{}
	MeshData(const MeshData& m) : mesh_(m.mesh_)
	{}

	MeshData(const Mesh* mesh) : mesh_(mesh)
	{}
	
	void draw(rendering::DrawingType primitive)
	{
		if (!render_.is_primitive_uptodate(primitive))
			render_.init_primitives(*mesh_, primitive);
		render_.draw(primitive);
	}

	void init_primitives(rendering::DrawingType primitive)
	{
		render_.init_primitives(*mesh_, primitive);
	}

	void update_bb(const Attribute<Vec3>* vertex_position)
	{
		for (uint32 i = 0; i < 3; ++i)
		{
			bb_min_[i] = std::numeric_limits<float64>::max();
			bb_max_[i] = std::numeric_limits<float64>::lowest();
		}
		for (const Vec3& v : *vertex_position)
		{
			for (uint32 i = 0; i < 3; ++i)
			{
				if (v[i] < bb_min_[i])
					bb_min_[i] = v[i];
				if (v[i] > bb_max_[i])
					bb_max_[i] = v[i];
			}
		}
	}

	rendering::VBO* vbo(const std::string& name)
	{
		if (auto it = vbos_.find(name); it != vbos_.end())
			return it->second.get();
		else
			return nullptr;
	}

	// rendering::VBO* create_vbo(const std::string& name)
	// {
	// 	rendering::VBO* v = vbo(name);
	// 	if (!v)
	// 	{
	// 		std::shared_ptr<Attribute<Vec3>> attribute3 = get_attribute<Vec3, Vertex>(*mesh_, name);
	// 		if (attribute3)
	// 		{
	// 			const auto [it, inserted] = vbos_.emplace(name, std::make_unique<rendering::VBO>(3));
	// 			v = it->second.get();
	// 			rendering::update_vbo<Vec3>(attribute3.get(), v);
	// 			return v;
	// 		}
	// 		std::shared_ptr<Attribute<Vec2>> attribute2 = get_attribute<Vec2, Vertex>(*mesh_, name);
	// 		if (attribute2)
	// 		{
	// 			const auto [it, inserted] = vbos_.emplace(name, std::make_unique<rendering::VBO>(2));
	// 			v = it->second.get();
	// 			rendering::update_vbo<Vec2>(attribute2.get(), v);
	// 			return v;
	// 		}
	// 		std::shared_ptr<Attribute<Scalar>> attribute1 = get_attribute<Scalar, Vertex>(*mesh_, name);
	// 		if (attribute1)
	// 		{
	// 			const auto [it, inserted] = vbos_.emplace(name, std::make_unique<rendering::VBO>(1));
	// 			v = it->second.get();
	// 			rendering::update_vbo<Scalar>(attribute1.get(), v);
	// 			return v;
	// 		}
	// 	}
	// 	return v;
	// }

	void update_vbo(const std::string& name)
	{
		rendering::VBO* v = vbo(name);

		std::shared_ptr<Attribute<Vec3>> attribute3 = get_attribute<Vec3, Vertex>(*mesh_, name);
		if (attribute3)
		{
			if (!v)
			{
				const auto [it, inserted] = vbos_.emplace(name, std::make_unique<rendering::VBO>(3));
				v = it->second.get();
			}
			rendering::update_vbo<Vec3>(attribute3.get(), v);
			return;
		}
		std::shared_ptr<Attribute<Vec2>> attribute2 = get_attribute<Vec2, Vertex>(*mesh_, name);
		if (attribute2)
		{
			if (!v)
			{
				const auto [it, inserted] = vbos_.emplace(name, std::make_unique<rendering::VBO>(2));
				v = it->second.get();
			}
			rendering::update_vbo<Vec2>(attribute2.get(), v);
			return;
		}
		std::shared_ptr<Attribute<Scalar>> attribute1 = get_attribute<Scalar, Vertex>(*mesh_, name);
		if (attribute1)
		{
			if (!v)
			{
				const auto [it, inserted] = vbos_.emplace(name, std::make_unique<rendering::VBO>(1));
				v = it->second.get();
			}
			rendering::update_vbo<Scalar>(attribute1.get(), v);
			return;
		}
	}

	Vec3 bb_min_, bb_max_;

private:

	const Mesh* mesh_;
	rendering::MeshRender render_;
	std::unordered_map<std::string, std::unique_ptr<rendering::VBO>> vbos_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_MESH_PROVIDER_MESH_DATA_H_

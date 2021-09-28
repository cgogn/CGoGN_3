/*******************************************************************************
 * CGoGN                                                                        *
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

#ifndef CGOGN_MODULE_SURFACE_FILTERING_H_
#define CGOGN_MODULE_SURFACE_FILTERING_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/geometry/algos/filtering.h>
#include <cgogn/geometry/algos/laplacian.h>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceFiltering : public Module
{
	static_assert(mesh_traits<MESH>::dimension >= 2, "SurfaceFiltering can only be used with meshes of dimension >= 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	SurfaceFiltering(const App& app)
		: Module(app, "SurfaceFiltering (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  selected_vertex_attribute_(nullptr)
	{
	}
	~SurfaceFiltering()
	{
	}

	void filter_mesh(MESH& m, Attribute<Vec3>* vertex_attribute)
	{
		std::shared_ptr<Attribute<Vec3>> filtered_vertex_attribute =
			add_attribute<Vec3, Vertex>(m, "__filtered_attribute");
		geometry::filter_average<Vec3>(m, vertex_attribute, filtered_vertex_attribute.get());
		vertex_attribute->swap(filtered_vertex_attribute.get());
		remove_attribute<Vertex>(m, filtered_vertex_attribute);

		mesh_provider_->emit_attribute_changed(m, vertex_attribute);
	}

	void regularize(MESH& m, Attribute<Vec3>* vertex_position)
	{
		auto vertex_index = add_attribute<uint32, Vertex>(m, "__vertex_index");
		auto vertex_position_laplacian = add_attribute<Vec3, Vertex>(m, "__vertex_position_laplacian");

		uint32 nb_vertices = 0;
		foreach_cell(m, [&](Vertex v) -> bool {
			value<uint32>(m, vertex_index, v) = nb_vertices++;
			return true;
		});

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> LAPL =
			geometry::cotan_laplacian_matrix(m, vertex_index.get(), vertex_position);
		Eigen::MatrixXd vpos(nb_vertices, 3);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			const Vec3& pv = value<Vec3>(m, vertex_position, v);
			uint32 vidx = value<uint32>(m, vertex_index, v);
			vpos(vidx, 0) = pv[0];
			vpos(vidx, 1) = pv[1];
			vpos(vidx, 2) = pv[2];
			return true;
		});
		Eigen::MatrixXd poslapl(nb_vertices, 3);
		poslapl = LAPL * vpos;
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			Vec3& vplapl = value<Vec3>(m, vertex_position_laplacian, v);
			uint32 vidx = value<uint32>(m, vertex_index, v);
			vplapl[0] = poslapl(vidx, 0);
			vplapl[1] = poslapl(vidx, 1);
			vplapl[2] = poslapl(vidx, 2);
			return true;
		});

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(2 * nb_vertices, nb_vertices);
		std::vector<Eigen::Triplet<Scalar>> Acoeffs;
		Acoeffs.reserve(nb_vertices * 10);
		foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			auto vertices = adjacent_vertices_through_edge(m, v);
			uint32 d = vertices.size();
			for (Vertex av : vertices)
			{
				uint32 avidx = value<uint32>(m, vertex_index, av);
				Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(avidx), 1));
			}
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx), int(vidx), -1 * Scalar(d)));
			Acoeffs.push_back(Eigen::Triplet<Scalar>(int(nb_vertices + vidx), int(vidx), 10));
			return true;
		});
		A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> solver(A);

		Eigen::MatrixXd x(nb_vertices, 3);
		Eigen::MatrixXd b(2 * nb_vertices, 3);

		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			const Vec3& l = value<Vec3>(m, vertex_position_laplacian, v);
			b(vidx, 0) = l[0];
			b(vidx, 1) = l[1];
			b(vidx, 2) = l[2];
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			b(nb_vertices + vidx, 0) = 10 * pos[0];
			b(nb_vertices + vidx, 1) = 10 * pos[1];
			b(nb_vertices + vidx, 2) = 10 * pos[2];
			x(vidx, 0) = pos[0];
			x(vidx, 1) = pos[1];
			x(vidx, 2) = pos[2];
			return true;
		});

		x = solver.solveWithGuess(b, x);

		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			Vec3& pos = value<Vec3>(m, vertex_position, v);
			pos[0] = x(vidx, 0);
			pos[1] = x(vidx, 1);
			pos[2] = x(vidx, 2);
			return true;
		});

		remove_attribute<Vertex>(m, vertex_index);
		remove_attribute<Vertex>(m, vertex_position_laplacian);

		mesh_provider_->emit_attribute_changed(m, vertex_position);
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
	}

	void interface() override
	{
		imgui_mesh_selector(mesh_provider_, selected_mesh_, "Surface", [&](MESH& m) {
			selected_mesh_ = &m;
			selected_vertex_attribute_.reset();
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_attribute_, "Attribute",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_attribute_ = attribute; });

			if (selected_vertex_attribute_)
			{
				if (ImGui::Button("Filter"))
					filter_mesh(*selected_mesh_, selected_vertex_attribute_.get());
				if (ImGui::Button("Regularize"))
					regularize(*selected_mesh_, selected_vertex_attribute_.get());
			}
		}
	}

private:
	MESH* selected_mesh_;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_attribute_;
	MeshProvider<MESH>* mesh_provider_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_FILTERING_H_

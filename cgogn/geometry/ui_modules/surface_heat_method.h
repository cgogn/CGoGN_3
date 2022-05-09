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

#ifndef CGOGN_MODULE_SURFACE_HEAT_METHOD_H_
#define CGOGN_MODULE_SURFACE_HEAT_METHOD_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <queue>
#include <unordered_set>
#include <utility>

namespace cgogn
{

namespace ui
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;

Scalar vertex_gradient_divergence(const CMap2& m, CMap2::Vertex v, const CMap2::Attribute<Vec3>* face_gradient,
								  const CMap2::Attribute<Vec3>* vertex_position)
{
	Scalar div = 0.0;
	std::vector<CMap2::Edge> edges = incident_edges(m, v);
	for (uint32 i = 0; i < edges.size(); ++i)
	{
		CMap2::Edge e1 = edges[i];
		CMap2::Edge e2 = edges[(i + 1) % edges.size()];

		CMap2::Face f(e1.dart);

		const Vec3& X = value<Vec3>(m, face_gradient, f);

		const Vec3& p0 = value<Vec3>(m, vertex_position, v);
		const Vec3& p1 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi1(m, e1.dart)));
		const Vec3& p2 = value<Vec3>(m, vertex_position, CMap2::Vertex(phi1(m, e2.dart)));

		Vec3 vecR = p0 - p2;
		Vec3 vecL = p1 - p2;
		Scalar cotValue1 = vecR.dot(vecL) / vecR.cross(vecL).norm();

		vecR = p2 - p1;
		vecL = p0 - p1;
		Scalar cotValue2 = vecR.dot(vecL) / vecR.cross(vecL).norm();

		div += cotValue1 * (p1 - p0).dot(X) + cotValue2 * (p2 - p0).dot(X);
	}
	return div / 2.0;
}

template <typename MESH>
class SurfaceHeatMethod : public Module
{
	static_assert(mesh_traits<MESH>::dimension == 2, "SurfaceHeatMethod can only be used with meshes of dimension 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;
	using Edge = typename mesh_traits<MESH>::Edge;

public:
	SurfaceHeatMethod(const App& app) : Module(app, "SurfaceHeatMethod (" + std::string{mesh_traits<MESH>::name} + ")")
	{
	}

	~SurfaceHeatMethod()
	{
	}

	void euclidean_distance(const MESH& m, const Attribute<Vec3>* vertex_position,
							const CellsSet<MESH, Vertex>* source_vertices, Attribute<Scalar>* vertex_euclidean_dist)
	{
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			Scalar s = std::numeric_limits<Scalar>::max();
			const Vec3& pos = value<Vec3>(m, vertex_position, v);
			source_vertices->foreach_cell([&](Vertex vec) {
				Scalar d = (pos - value<Vec3>(m, vertex_position, vec)).norm();
				if (s > d)
					s = d;
			});
			value<Scalar>(m, vertex_euclidean_dist, v) = s;
			return true;
		});

		mesh_provider_->emit_attribute_changed(m, vertex_euclidean_dist);
	}

	void geodesic_distance(MESH& m, const Attribute<Vec3>* vertex_position,
						   const CellsSet<MESH, Vertex>* source_vertices, Attribute<Scalar>* vertex_geodesic_distance)
	{
		auto vertex_index = add_attribute<uint32, Vertex>(m, "__vertex_index");

		uint32 nb_vertices = 0;
		foreach_cell(m, [&](Vertex v) -> bool {
			value<uint32>(m, vertex_index, v) = nb_vertices++;
			return true;
		});

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Lc =
			geometry::cotan_operator_matrix(m, vertex_index.get(), vertex_position);

		auto vertex_area = add_attribute<Scalar, Vertex>(m, "__vertex_area");
		geometry::compute_area<Vertex>(m, vertex_position, vertex_area.get());

		Eigen::VectorXd A(nb_vertices);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			A(vidx) = value<Scalar>(m, vertex_area, v);
			return true;
		});

		Eigen::VectorXd u0(nb_vertices);
		u0.setZero();
		source_vertices->foreach_cell([&](Vertex v) {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			u0(vidx) = 1.0;
		});

		Scalar h = geometry::mean_edge_length(m, vertex_position);
		Scalar t = h * h;

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Am(A.asDiagonal());
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> heat_solver(Am - t * Lc);
		Eigen::VectorXd u = heat_solver.solve(u0);

		auto vertex_heat = get_or_add_attribute<Scalar, Vertex>(m, "__vertex_heat");
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			value<Scalar>(m, vertex_heat, v) = u(vidx);
			return true;
		});

		auto face_heat_gradient = get_or_add_attribute<Vec3, Face>(m, "__face_heat_gradient");
		parallel_foreach_cell(m, [&](Face f) -> bool {
			Vec3 g(0, 0, 0);
			Vec3 n = geometry::normal(m, f, vertex_position);
			Scalar a = geometry::area(m, f, vertex_position);
			std::vector<Vertex> vertices = incident_vertices(m, f);
			for (uint32 i = 0; i < vertices.size(); ++i)
			{
				Vec3 e = value<Vec3>(m, vertex_position, vertices[(i + 2) % vertices.size()]) -
						 value<Vec3>(m, vertex_position, vertices[(i + 1) % vertices.size()]);
				g += value<Scalar>(m, vertex_heat, vertices[i]) * n.cross(e);
			}
			g /= 2 * a;
			value<Vec3>(m, face_heat_gradient, f) = -1.0 * g.normalized();
			return true;
		});

		auto vertex_heat_gradient_div = get_or_add_attribute<Scalar, Vertex>(m, "__vertex_heat_gradient_div");
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			Scalar d = vertex_gradient_divergence(m, v, face_heat_gradient.get(), vertex_position);
			value<Scalar>(m, vertex_heat_gradient_div, v) = d;
			return true;
		});

		Eigen::VectorXd b(nb_vertices);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			b(vidx) = value<Scalar>(m, vertex_heat_gradient_div, v);
			return true;
		});

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> poisson_solver(Lc);
		Eigen::VectorXd dist = poisson_solver.solve(b);

		Scalar min = dist.minCoeff();
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			value<Scalar>(m, vertex_geodesic_distance, v) = dist(vidx) - min;
			return true;
		});

		remove_attribute<Vertex>(m, vertex_index);
		remove_attribute<Vertex>(m, vertex_area);
		// remove_attribute<Vertex>(m, vertex_heat);
		// remove_attribute<Face>(m, face_heat_gradient);
		// remove_attribute<Vertex>(m, vertex_heat_gradient_div);

		mesh_provider_->emit_attribute_changed(m, vertex_geodesic_distance);
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
			selected_vertex_position_.reset();
			mesh_provider_->mesh_data(m).outlined_until_ = App::frame_time_ + 1.0;
		});

		if (selected_mesh_)
		{
			MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);

			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_position_, "Vertex Position",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position_ = attribute; });

			imgui_combo_cells_set(md, selected_vertices_set_, "Source vertices",
								  [&](CellsSet<MESH, Vertex>* cs) { selected_vertices_set_ = cs; });

			if (selected_vertex_position_ && selected_vertices_set_)
			{
				imgui_combo_attribute<Vertex, Scalar>(*selected_mesh_, euclidean_distance_vertex_,
													  "Vertex Euclidean distance",
													  [&](const std::shared_ptr<Attribute<Scalar>>& attribute) {
														  euclidean_distance_vertex_ = attribute;
													  });

				if (ImGui::Button("Compute Euclidean Distance"))
				{
					if (!euclidean_distance_vertex_)
						euclidean_distance_vertex_ =
							get_or_add_attribute<Scalar, Vertex>(*selected_mesh_, "euclidean distance");
					euclidean_distance(*selected_mesh_, selected_vertex_position_.get(), selected_vertices_set_,
									   euclidean_distance_vertex_.get());
				}

				imgui_combo_attribute<Vertex, Scalar>(*selected_mesh_, geodesic_distance_vertex_,
													  "Vertex Geodesic distance",
													  [&](const std::shared_ptr<Attribute<Scalar>>& attribute) {
														  geodesic_distance_vertex_ = attribute;
													  });
				if (ImGui::Button("Compute Geodesic Distance"))
				{
					if (!geodesic_distance_vertex_)
						geodesic_distance_vertex_ =
							get_or_add_attribute<Scalar, Vertex>(*selected_mesh_, "geodesic distance");
					geodesic_distance(*selected_mesh_, selected_vertex_position_.get(), selected_vertices_set_,
									  geodesic_distance_vertex_.get());
				}

				static double t_multiplier = 1.0;
				ImGui::InputDouble("Scalar t_multiplier", &t_multiplier, 0.01f, 100.0f, "%.3f");
			}
		}
	}

private:
	MESH* selected_mesh_ = nullptr;

	// geometry::DistanceHeatSolver<MESH>* distanceHeatSolver = nullptr;

	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_ = nullptr;
	CellsSet<MESH, Vertex>* selected_vertices_set_ = nullptr;

	MeshProvider<MESH>* mesh_provider_ = nullptr;

	std::shared_ptr<Attribute<Scalar>> euclidean_distance_vertex_ = nullptr;
	std::shared_ptr<Attribute<Scalar>> geodesic_distance_vertex_ = nullptr;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_HEAT_METHOD_H_

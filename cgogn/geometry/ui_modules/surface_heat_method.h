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

#include <chrono>
#include <utility>
#include <bits/stdc++.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <limits>
#include <cmath>
#include <queue>
#include <iostream>
#include <unordered_set>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <Eigen/Dense>

#include <Eigen/Sparse>

#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/algos/laplacian.h>
#include <cgogn/core/functions/attributes.h>


#include <cgogn/geometry/algos/vector_heat_solver.h>
#include <cgogn/geometry/algos/distance_heat_solver.h>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceHeatMethod : public Module
{
	static_assert(mesh_traits<MESH>::dimension == 2, "SurfaceHeatMethod can only be used with meshes of dimension 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;
	using Edge = typename mesh_traits<MESH>::Edge;

	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

public:
	SurfaceHeatMethod(const App& app)
		: Module(app, "SurfaceHeatMethod (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  selected_vertex_position_(nullptr), selected_vertices_set_(nullptr), euclidean_dist_vert(nullptr), euclidean_dist_face(nullptr),
		  geo_dist_vert(nullptr), geo_dist_face(nullptr), heat_distance(nullptr), selected_face_normal_(nullptr)
	{
	}
	~SurfaceHeatMethod()
	{
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
	}

	void euclid_dist(){
		if (euclidean_dist_vert == nullptr) {
			euclidean_dist_vert = get_or_add_attribute<Scalar,Vertex>(*selected_mesh_,"euclidean_dist_vert");
		}
		if (euclidean_dist_face == nullptr) {
			euclidean_dist_face = get_or_add_attribute<Scalar,Face>(*selected_mesh_,"euclidean_dist_face");
		}
		foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
			double s = std::numeric_limits<double>::max();
			const Vec3& pos = value<Vec3>(*selected_mesh_, selected_vertex_position_, v);
			
			selected_vertices_set_->foreach_cell([&, pos](Vertex vec) {
				double d = (pos - value<Vec3>(*selected_mesh_, selected_vertex_position_, vec)).norm();
				if (s > d){
					s = d;
				}
			});
			value<Scalar>(*selected_mesh_, euclidean_dist_vert, v) = s;
			
			//value<Scalar>(*selected_mesh_, euclidean_dist_vert, v) = pos.norm();
			return true;
		});
		foreach_cell(*selected_mesh_, [&](Face f) -> bool {
			double v = 0;
			std::vector<Vertex> vertices = incident_vertices(*selected_mesh_, f);
			for (int i = 0; i < vertices.size(); i++){
				v += value<Scalar>(*selected_mesh_, euclidean_dist_vert, vertices[i]);
			}
			/*
			v = value<Scalar>(*selected_mesh_, euclidean_dist_vert, vertices[0])
					+ value<Scalar>(*selected_mesh_, euclidean_dist_vert, vertices[1])
					+ value<Scalar>(*selected_mesh_, euclidean_dist_vert, vertices[2]);
					*/
			value<Scalar>(*selected_mesh_, euclidean_dist_face, f) = v/vertices.size();

			return true;
		});
	}

	void geo_dist(){

		if (geo_dist_vert == nullptr) {
			geo_dist_vert = get_or_add_attribute<Scalar,Vertex>(*selected_mesh_,"geo_dist_vert");
		}
		if (geo_dist_face == nullptr) {
			geo_dist_face = get_or_add_attribute<Scalar,Face>(*selected_mesh_,"geo_dist_face");
		}

		double max_value = std::numeric_limits<double>::max();
		foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
			value<Scalar>(*selected_mesh_, geo_dist_vert, v) = max_value;
			return true;
		});
		std::queue<Vertex> tmp;


		//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
		selected_vertices_set_->foreach_cell([&](Vertex vec) {
			value<Scalar>(*selected_mesh_, geo_dist_vert, vec) = 0.0;
			tmp.push(vec);
		});

		while (!tmp.empty()){
			Vertex v = tmp.front();
			tmp.pop();
			const Vec3& pos = value<Vec3>(*selected_mesh_, selected_vertex_position_, v);
			double point_value = value<Scalar>(*selected_mesh_, geo_dist_vert, v);
			foreach_adjacent_vertex_through_edge(*selected_mesh_, v, [&, pos, point_value](Vertex vert) -> bool{
				double d = point_value + (pos - value<Vec3>(*selected_mesh_, selected_vertex_position_, vert)).norm();
				if (value<Scalar>(*selected_mesh_, geo_dist_vert, vert) > d){
					value<Scalar>(*selected_mesh_, geo_dist_vert, vert) = d;
						tmp.push(vert);
				}
				return true;
			});
		}

		foreach_cell(*selected_mesh_, [&](Face f) -> bool {
			double v = 0;
			std::vector<Vertex> vertices = incident_vertices(*selected_mesh_, f);
			for (int i = 0; i < vertices.size(); i++){
				v += value<Scalar>(*selected_mesh_, geo_dist_vert, vertices[i]);
			}
			value<Scalar>(*selected_mesh_, geo_dist_face, f) = v/vertices.size();

			return true;
		});
	}

/*
	struct heat_precompute_data {
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>>* heat_solver;
	}
*/

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>>* compute_heat_solver(
				MESH& mesh,
				std::shared_ptr<Attribute<Vec3>> position,
				std::shared_ptr<Attribute<uint32>> vertex_index,
				uint32 nb_vertices,
				Eigen::SparseMatrix<Scalar, Eigen::ColMajor>* Lc,
				//std::shared_ptr<Attribute<Scalar>> area = nullptr,
				double time_multiplier = 1.0,
				std::vector<std::pair<std::string,std::chrono::steady_clock::time_point>>* time_logger = nullptr,
				bool keep_local_attribute = false)
	{

		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("start", std::chrono::steady_clock::now()));


		auto area = get_or_add_attribute<Scalar, Vertex>(mesh, "__area");
		geometry::compute_area<Vertex>(mesh, position.get(), area.get());

		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("init", std::chrono::steady_clock::now()));

		*Lc = geometry::cotan_operator_matrix(mesh, vertex_index.get(), position.get());

		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("compute Lc", std::chrono::steady_clock::now()));

		Eigen::VectorXd areas = Eigen::VectorXd(nb_vertices);


		foreach_cell(mesh, [&](Vertex v) -> bool {
			areas(value<uint32>(mesh, vertex_index, v)) =
				value<Scalar>(mesh, area, v);
			return true;
		});

		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("construct areas matrix", std::chrono::steady_clock::now()));
		double t = 0;
		int nb_edges = 0;
		foreach_cell(mesh, [&](Edge e) -> bool {
			t += geometry::length(mesh, e, position.get());
			nb_edges++;
			return true;
		});

		t /= nb_edges;
		t *= t;
		t *= time_multiplier;

		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("compute t", std::chrono::steady_clock::now()));

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>>* solver = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>>(Eigen::SparseMatrix<Scalar, Eigen::ColMajor>(areas.asDiagonal()) - (t*(*Lc)));
		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("construct heat solver", std::chrono::steady_clock::now()));
		if (!keep_local_attribute){
			remove_attribute<Vertex>(mesh, area);
		}
		return solver;
	}

	void heat_compute(MESH& mesh,
				std::shared_ptr<Attribute<Scalar>> out_heat_distance,
				std::shared_ptr<Attribute<Vec3>> position,
				std::shared_ptr<Attribute<Vec3>> face_normal,
				std::shared_ptr<Attribute<uint32>> vertex_index,
				uint32 nb_vertices,
				Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>>& solver,
				CellsSet<MESH, Vertex>& cell_source,
				Eigen::SparseMatrix<Scalar, Eigen::ColMajor>& Lc,
				std::vector<std::pair<std::string,std::chrono::steady_clock::time_point>>* time_logger = nullptr,
				bool keep_local_attribute = false){

		Eigen::VectorXd source = Eigen::VectorXd::Zero(nb_vertices);

		cell_source.foreach_cell([&](Vertex vec) {
			source(value<uint32>(mesh, vertex_index.get(), vec)) = 1;
		});


		auto u = solver.solve(source);
		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("solve heat diffusion", std::chrono::steady_clock::now()));

		auto heat_attribute = get_or_add_attribute<Scalar,Vertex>(mesh,"__heat");
		foreach_cell(mesh, [&](Vertex v) -> bool {
			value<Scalar>(mesh, heat_attribute, v) = u(value<uint32>(mesh, vertex_index, v));
			return true;
		});

		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("apply heat to mesh", std::chrono::steady_clock::now()));
				
		// gradient
		auto gradient_attribute = get_or_add_attribute<Vec3, Face>(mesh, "__heat_gradient");
		auto gradient_normalized_attribute = get_or_add_attribute<Vec3, Face>(mesh, "__heat_normalized_gradient");
		foreach_cell(mesh, [&](Face f) -> bool {
			Vec3 Ni = value<Vec3>(mesh, face_normal, f);
			std::vector<Vertex> vertices = incident_vertices(mesh, f);
			Vec3 grad;
			for (int i = 0; i < 3; i++){
				Scalar Ui = value<Scalar>(mesh, heat_attribute, vertices[(i+2)%3]);
				Vec3 ei = value<Vec3>(mesh, position, vertices[(i+1)%3])
						- value<Vec3>(mesh, position, vertices[i]);
				grad += Ni.cross(ei) * Ui;
			}
			grad /= 2.0 * geometry::area(mesh, f, position.get());
			value<Vec3>(mesh, gradient_attribute, f) = grad;
			value<Vec3>(mesh, gradient_normalized_attribute, f) = (grad / grad.norm() * -1.0);
			return true;
		});
		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("compute gradient", std::chrono::steady_clock::now()));

		auto integrated_divergence = get_or_add_attribute<Scalar,Vertex>(mesh, "__integrated_divergence");
		foreach_cell(mesh, [&](Vertex v) -> bool {
			Vec3 vpos = value<Vec3>(mesh, position, v);
			Scalar div = 0.0;
			foreach_incident_face(mesh, v, [&](Face f) -> bool {
				Vec3 Xj = value<Vec3>(mesh, gradient_normalized_attribute, f);
				std::vector<Vertex> vertices = incident_vertices(mesh, f);
				std::remove_if(vertices.begin(), vertices.end(), [&](Vertex item) -> bool {
					return value<uint32_t>(mesh, vertex_index, item) == value<uint32_t>(mesh, vertex_index, v);
				});
				Vec3 p0 = value<Vec3>(mesh, position, vertices[0]);
				Vec3 p1 = value<Vec3>(mesh, position, vertices[1]);
				Vec3 e1 = p0 - vpos;
				Vec3 e2 = p1 - vpos;
				double teta1 = geometry::angle(-e2, p0 - p1);
				double teta2 = geometry::angle(-e1, p1 - p0);
				div += 1.0/tan(teta1) * e1.dot(Xj) + 1.0/tan(teta2) * e2.dot(Xj);

				return true;
			});
			value<Scalar>(mesh, integrated_divergence, v) = div * 0.5;
			return true;
		});

		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("compute integrated divergence", std::chrono::steady_clock::now()));

		Eigen::VectorXd b = Eigen::VectorXd::Zero(nb_vertices);
		foreach_cell(mesh, [&](Vertex v) -> bool {
			b(value<uint32>(mesh, vertex_index.get(), v)) = 
				value<Scalar>(mesh, integrated_divergence, v);
			return true;
		});

		auto geo_solver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>>(Eigen::SparseMatrix<Scalar, Eigen::ColMajor>(Lc));

		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("construct geodesic solver", std::chrono::steady_clock::now()));

		auto phi = geo_solver.solve(b);

		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("solve geodesic distance", std::chrono::steady_clock::now()));

		
		//out_heat_distance = get_or_add_attribute<Scalar, Vertex>(mesh, "heat_distance");
		parallel_foreach_cell(mesh, [&](Vertex v) -> bool {
			value<Scalar>(mesh, out_heat_distance, v) =
				phi(value<uint32>(mesh, vertex_index.get(), v));
			return true;
		});

		Scalar min = std::numeric_limits<Scalar>::max();
		foreach_cell(mesh, [&](Vertex v) -> bool {
			if (value<Scalar>(mesh, out_heat_distance, v) < min){
				min = value<Scalar>(mesh, out_heat_distance, v);
			}
			return true;
		});
		parallel_foreach_cell(mesh, [&](Vertex v) -> bool {
			value<Scalar>(mesh, out_heat_distance, v) -= min;
			return true;
		});

		if (time_logger != nullptr) time_logger->push_back(std::pair<std::string,std::chrono::steady_clock::time_point>
				("apply geodesic distance", std::chrono::steady_clock::now()));

/*
		mesh_provider_->emit_attribute_changed(mesh, gradient_attribute.get());
		mesh_provider_->emit_attribute_changed(mesh, gradient_normalized_attribute.get());
		mesh_provider_->emit_attribute_changed(mesh, heat_distance.get());
		mesh_provider_->emit_attribute_changed(mesh, integrated_divergence.get());
		*/
		if (!keep_local_attribute){
			remove_attribute<Vertex>(mesh, gradient_attribute);
			remove_attribute<Vertex>(mesh, gradient_normalized_attribute);
			remove_attribute<Vertex>(mesh, integrated_divergence);
			remove_attribute<Vertex>(mesh, heat_attribute);
		}
	}

	void heat_compute_one_step(MESH& mesh,
				std::shared_ptr<Attribute<Scalar>> out_heat_distance,
				std::shared_ptr<Attribute<Vec3>> position,
				std::shared_ptr<Attribute<Vec3>> face_normal,
				CellsSet<MESH, Vertex>& cell_source,
				double time_multiplier = 1.0,
				std::vector<std::pair<std::string,std::chrono::steady_clock::time_point>>* time_logger = nullptr,
				bool keep_local_attribute = false
				)
	{
		auto vertex_index = get_or_add_attribute<uint32, Vertex>(mesh, "__vertex_index");

		uint32 nb_vertices = 0;
		foreach_cell(mesh, [&](Vertex v) -> bool {
			value<uint32>(mesh, vertex_index, v) = nb_vertices++;
			return true;
		});
		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Lc;
		auto solver = compute_heat_solver(mesh, position, vertex_index, nb_vertices, &Lc, time_multiplier, time_logger, keep_local_attribute);

		heat_compute(mesh, out_heat_distance, position, face_normal, vertex_index, nb_vertices, *solver, cell_source, Lc, time_logger, keep_local_attribute);
		//heat_compute(mesh, out_heat_distance, position, face_normal, vertex_index, nb_vertices, *solver, cell_source, Lc, time_multiplier, time_logger, keep_local_attribute);
		if (!keep_local_attribute){
			remove_attribute<Vertex>(mesh, vertex_index);
		}
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
			heat_distance = get_or_add_attribute<Scalar,Vertex>(*selected_mesh_,"heat_distance");

			imgui_combo_attribute<Vertex, Vec3>(
				*selected_mesh_, selected_vertex_position_, "Vertex Position",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_vertex_position_ = attribute; });

			imgui_combo_attribute<Face, Vec3>(
				*selected_mesh_, selected_face_normal_, "Face Normal",
				[&](const std::shared_ptr<Attribute<Vec3>>& attribute) { selected_face_normal_ = attribute; });
/*
			imgui_combo_attribute<Vertex, uint32>(
				*selected_mesh_, selected_vertex_index_, "Vertex Index",
				[&](const std::shared_ptr<Attribute<uint32>>& attribute) { selected_vertex_index_ = attribute; });
				*/

			imgui_combo_cells_set(md, selected_vertices_set_, "Source vertices",
								  [&](CellsSet<MESH, Vertex>* cs) { selected_vertices_set_ = cs; });
/*
			static char index_name [50] = "vertex_index";
			static uint32 nb_vertices = 0;
			ImGui::InputText("Index Name", index_name, sizeof(index_name));
			if (ImGui::Button("Generate Index")){
				auto vertex_index = get_or_add_attribute<uint32, Vertex>(*selected_mesh_, index_name);

				nb_vertices = 0;
				foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
					value<uint32>(*selected_mesh_, vertex_index, v) = nb_vertices++;
					return true;
				});
			}
			*/
			if (selected_vertex_position_ != nullptr && selected_vertices_set_ != nullptr)
			{
				if (ImGui::Button("Compute Euclidean Distance")){
					euclid_dist();
					mesh_provider_->emit_attribute_changed(*selected_mesh_, euclidean_dist_face.get());
				}
				if (ImGui::Button("Compute Topo Distance")){
					geo_dist();
					mesh_provider_->emit_attribute_changed(*selected_mesh_, geo_dist_face.get());
				}
				static double t_multiplier = 1.0;
				ImGui::InputDouble("Scalar t_multiplier",
											  &t_multiplier, 0.01f, 100.0f,
											  "%.3f");
				static std::vector<std::pair<std::string,std::chrono::steady_clock::time_point>> time_log;
				if (selected_face_normal_ != nullptr){
					{
						if (ImGui::Button("Compute Heat With Timelog")){
							time_log.clear();
							heat_compute_one_step(*selected_mesh_, heat_distance, selected_vertex_position_, selected_face_normal_, *selected_vertices_set_, t_multiplier, &time_log, true);
							mesh_provider_->emit_attribute_changed(*selected_mesh_, heat_distance.get());
						}
						if (ImGui::Button("(Re)Generate DistanceHeatSolver")){
							if (this->distanceHeatSolver != nullptr){
								delete distanceHeatSolver;
							}
							time_log.clear();
							time_log.push_back(std::pair<std::string,std::chrono::steady_clock::time_point>("start", std::chrono::steady_clock::now()));
							distanceHeatSolver = new geometry::DistanceHeatSolver<MESH>(*selected_mesh_, selected_vertex_position_, selected_face_normal_, t_multiplier);
							time_log.push_back(std::pair<std::string,std::chrono::steady_clock::time_point>("(Re)Generate DistanceHeatSolver", std::chrono::steady_clock::now()));
						}
						if (distanceHeatSolver != nullptr){
							if (ImGui::Button("Compute From Source")){
								if (heat_distance == nullptr){
									heat_distance = get_or_add_attribute<Scalar, Vertex>(*selected_mesh_, "heat_distance");
								}
								time_log.clear();
								time_log.push_back(std::pair<std::string,std::chrono::steady_clock::time_point>("start", std::chrono::steady_clock::now()));
								distanceHeatSolver->solve(heat_distance.get(), *selected_vertices_set_);
								time_log.push_back(std::pair<std::string,std::chrono::steady_clock::time_point>("Compute From Source", std::chrono::steady_clock::now()));
								mesh_provider_->emit_attribute_changed(*selected_mesh_, heat_distance.get());
							}
							if (ImGui::Button("delete DistanceHeatSolver")){
								delete distanceHeatSolver;
								distanceHeatSolver = nullptr;
							}
						}
						/*
						if (ImGui::Button("Compute heat one step")){
							time_log.clear();
							heat_compute_one_step(*selected_mesh_, heat_distance, selected_vertex_position_, selected_face_normal_, *selected_vertices_set_, t_multiplier, &time_log, true);
							mesh_provider_->emit_attribute_changed(*selected_mesh_, heat_distance.get());
						}
						if (selected_vertex_index_ != nullptr && nb_vertices != 0){
							if (ImGui::Button("Compute Heat")){
								time_log.clear();
								Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Lc;
								auto solver = compute_heat_solver(*selected_mesh_, selected_vertex_position_, selected_vertex_index_, nb_vertices, &Lc, t_multiplier, &time_log, false);

								heat_compute(*selected_mesh_, heat_distance, selected_vertex_position_, selected_face_normal_,selected_vertex_index_, nb_vertices, *solver, *selected_vertices_set_, Lc, &time_log, false);
								mesh_provider_->emit_attribute_changed(*selected_mesh_, heat_distance.get());
								delete solver;
								geometry::VectorHeatSolver(*selected_mesh_, selected_vertex_position_.get());
								//geometry::VectorHeatSolver()
							}
						}
						*/
					}
				}
				
				if (time_log.size() > 0){
					if (ImGui::BeginTable("heat_time_log", 2, ImGuiTableFlags_Resizable)){
						for (int i = 1; i < time_log.size(); i++){
							ImGui::TableNextRow();
							ImGui::TableNextColumn();
							ImGui::Text(time_log[i].first.c_str());
							ImGui::TableNextColumn();
							std::string dtext = std::to_string(
								std::chrono::duration_cast<std::chrono::milliseconds>(time_log[i].second - time_log[i-1].second).count()
								)+ std::string(" ms");
							ImGui::Text(dtext.c_str());     
						}           
						ImGui::EndTable();
					}
				}
				
			}
		}
	}

private:
	MESH* selected_mesh_;
	geometry::DistanceHeatSolver<MESH>* distanceHeatSolver = nullptr;
	std::shared_ptr<Attribute<Vec3>> selected_vertex_position_;
	std::shared_ptr<Attribute<Vec3>> selected_face_normal_;
	CellsSet<MESH, Vertex>* selected_vertices_set_;
	MeshProvider<MESH>* mesh_provider_;
	std::shared_ptr<Attribute<Scalar>> euclidean_dist_vert;
	std::shared_ptr<Attribute<Scalar>> euclidean_dist_face;
	std::shared_ptr<Attribute<Scalar>> geo_dist_vert;
	std::shared_ptr<Attribute<Scalar>> geo_dist_face;
	std::shared_ptr<Attribute<Scalar>> heat_distance;
	std::shared_ptr<Attribute<uint32>> selected_vertex_index_ = nullptr;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_HEAT_METHOD_H_


/*

SimplicialDLT

filtering.h:71
laplacian.h:cotan_operator_matrix
edges : lengths.h
area.h : 
*/
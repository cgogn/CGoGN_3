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

#ifndef CGOGN_MODULE_SURFACE_DEFORMATION_H_
#define CGOGN_MODULE_SURFACE_DEFORMATION_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/modules/mesh_provider/mesh_provider.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/angle.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <GLFW/glfw3.h>

#include <Eigen/Sparse>
#include <boost/synapse/connect.hpp>
#include <memory>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class SurfaceDeformation : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension == 2, "SurfaceDeformation can only be used with meshes of dimension 2");

	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;

	using Vec3 = geometry::Vec3;
	using Mat3 = geometry::Mat3;
	using Scalar = geometry::Scalar;

	struct Parameters
	{
		Parameters()
			: vertex_position_(nullptr), selected_free_vertices_set_(nullptr), selected_handle_vertices_set_(nullptr),
			  initialized_(false), solver_ready_(false), vertex_position_init_(nullptr), vertex_diff_coord_(nullptr),
			  vertex_bi_diff_coord_(nullptr), vertex_rotation_matrix_(nullptr), vertex_rotated_diff_coord_(nullptr),
			  vertex_rotated_bi_diff_coord_(nullptr), vertex_index_(nullptr), edge_weight_(nullptr), solver_(nullptr)
		{
		}

		~Parameters()
		{
			if (solver_)
				delete solver_;
		}

		CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

		std::shared_ptr<Attribute<Vec3>> vertex_position_;

		CellsSet<MESH, Vertex>* selected_free_vertices_set_;
		CellsSet<MESH, Vertex>* selected_handle_vertices_set_;

		std::shared_ptr<boost::synapse::connection> cells_set_connection_;

		bool initialized_;
		bool solver_ready_;

		std::shared_ptr<Attribute<Vec3>> vertex_position_init_;
		std::shared_ptr<Attribute<Vec3>> vertex_diff_coord_;
		std::shared_ptr<Attribute<Vec3>> vertex_bi_diff_coord_;
		std::shared_ptr<Attribute<Mat3>> vertex_rotation_matrix_;
		std::shared_ptr<Attribute<Vec3>> vertex_rotated_diff_coord_;
		std::shared_ptr<Attribute<Vec3>> vertex_rotated_bi_diff_coord_;
		std::shared_ptr<Attribute<uint32>> vertex_index_;

		std::shared_ptr<Attribute<Scalar>> edge_weight_;

		std::unique_ptr<CellCache<MESH>> working_cells_;
		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> working_LAPL_;
		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> working_BILAPL_;

		Eigen::SparseLU<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>>* solver_;
	};

public:
	SurfaceDeformation(const App& app)
		: ViewModule(app, "SurfaceDeformation (" + std::string{mesh_traits<MESH>::name} + ")"), selected_mesh_(nullptr),
		  dragging_(false)
	{
	}
	~SurfaceDeformation()
	{
	}

private:
	void init_mesh(MESH* m)
	{
		Parameters& p = parameters_[m];
		p.working_cells_ = std::make_unique<CellCache<MESH>>(*m);
		p.cells_set_connection_ =
			boost::synapse::connect<typename MeshProvider<MESH>::template cells_set_changed<Vertex>>(
				m, [this, m](CellsSet<MESH, Vertex>* set) {
					Parameters& p = parameters_[m];
					if (p.selected_free_vertices_set_ == set || p.selected_handle_vertices_set_ == set)
						p.solver_ready_ = false;
				});
	}

	void initialize_mesh_data(MESH* m)
	{
		Parameters& p = parameters_[m];

		p.vertex_position_init_ = get_attribute<Vec3, Vertex>(*m, "position_init");
		if (!p.vertex_position_init_)
			p.vertex_position_init_ = add_attribute<Vec3, Vertex>(*m, "position_init");

		p.vertex_diff_coord_ = get_attribute<Vec3, Vertex>(*m, "diff_coord");
		if (!p.vertex_diff_coord_)
			p.vertex_diff_coord_ = add_attribute<Vec3, Vertex>(*m, "diff_coord");

		p.vertex_bi_diff_coord_ = get_attribute<Vec3, Vertex>(*m, "bi_diff_coord");
		if (!p.vertex_bi_diff_coord_)
			p.vertex_bi_diff_coord_ = add_attribute<Vec3, Vertex>(*m, "bi_diff_coord");

		p.vertex_rotation_matrix_ = get_attribute<Mat3, Vertex>(*m, "vertex_rotation_matrix");
		if (!p.vertex_rotation_matrix_)
			p.vertex_rotation_matrix_ = add_attribute<Mat3, Vertex>(*m, "vertex_rotation_matrix");

		p.vertex_rotated_diff_coord_ = get_attribute<Vec3, Vertex>(*m, "rotated_diff_coord");
		if (!p.vertex_rotated_diff_coord_)
			p.vertex_rotated_diff_coord_ = add_attribute<Vec3, Vertex>(*m, "rotated_diff_coord");

		p.vertex_rotated_bi_diff_coord_ = get_attribute<Vec3, Vertex>(*m, "rotated_bi_diff_coord");
		if (!p.vertex_rotated_bi_diff_coord_)
			p.vertex_rotated_bi_diff_coord_ = add_attribute<Vec3, Vertex>(*m, "rotated_bi_diff_coord");

		p.edge_weight_ = get_attribute<Scalar, Edge>(*m, "edge_weight");
		if (!p.edge_weight_)
			p.edge_weight_ = add_attribute<Scalar, Edge>(*m, "edge_weight");

		p.vertex_index_ = get_attribute<uint32, Vertex>(*m, "vertex_index");
		if (!p.vertex_index_)
			p.vertex_index_ = add_attribute<uint32, Vertex>(*m, "vertex_index");

		// initialize position init values
		p.vertex_position_init_->copy(p.vertex_position_.get());

		// initialize vertex rotation matrix
		Mat3 rm;
		rm.setZero();
		p.vertex_rotation_matrix_->fill(rm);

		// compute edges weight
		parallel_foreach_cell(*m, [&](Edge e) -> bool {
			std::vector<Scalar> angles = geometry::opposite_angles(*m, e, p.vertex_position_.get());
			Scalar& weight = value<Scalar>(*m, p.edge_weight_, e);
			for (Scalar a : angles)
				weight += std::tan(M_PI_2 - a);
			weight /= uint32(angles.size());
			return true;
		});

		// compute vertices laplacian
		uint32 nb_vertices = 0;
		foreach_cell(*m, [&](Vertex v) -> bool {
			value<uint32>(*m, p.vertex_index_, v) = nb_vertices++;
			return true;
		});

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> LAPL(nb_vertices, nb_vertices);
		std::vector<Eigen::Triplet<Scalar>> LAPLcoeffs;
		LAPLcoeffs.reserve(nb_vertices * 10);
		foreach_cell(*m, [&](Edge e) -> bool {
			Scalar w = value<Scalar>(*m, p.edge_weight_, e);
			auto vertices = incident_vertices(*m, e);
			uint32 vidx1 = value<uint32>(*m, p.vertex_index_, vertices[0]);
			uint32 vidx2 = value<uint32>(*m, p.vertex_index_, vertices[1]);
			LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx2), w));
			LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx1), w));
			LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx1), -w));
			LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx2), -w));
			return true;
		});
		LAPL.setFromTriplets(LAPLcoeffs.begin(), LAPLcoeffs.end());
		Eigen::MatrixXd vpos(nb_vertices, 3);
		parallel_foreach_cell(*m, [&](Vertex v) -> bool {
			const Vec3& pv = value<Vec3>(*m, p.vertex_position_, v);
			uint32 vidx = value<uint32>(*m, p.vertex_index_, v);
			vpos(vidx, 0) = pv[0];
			vpos(vidx, 1) = pv[1];
			vpos(vidx, 2) = pv[2];
			return true;
		});
		Eigen::MatrixXd lapl(nb_vertices, 3);
		lapl = LAPL * vpos;
		parallel_foreach_cell(*m, [&](Vertex v) -> bool {
			Vec3& dcv = value<Vec3>(*m, p.vertex_diff_coord_, v);
			uint32 vidx = value<uint32>(*m, p.vertex_index_, v);
			dcv[0] = lapl(vidx, 0);
			dcv[1] = lapl(vidx, 1);
			dcv[2] = lapl(vidx, 2);
			return true;
		});

		// compute vertices bi-laplacian
		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> BILAPL(nb_vertices, nb_vertices);
		BILAPL = LAPL * LAPL;
		Eigen::MatrixXd bilapl(nb_vertices, 3);
		bilapl = BILAPL * vpos;
		parallel_foreach_cell(*m, [&](Vertex v) -> bool {
			Vec3& bdcv = value<Vec3>(*m, p.vertex_bi_diff_coord_, v);
			uint32 vidx = value<uint32>(*m, p.vertex_index_, v);
			bdcv[0] = bilapl(vidx, 0);
			bdcv[1] = bilapl(vidx, 1);
			bdcv[2] = bilapl(vidx, 2);
			return true;
		});

		p.initialized_ = true;
		p.solver_ready_ = false;
	}

	void build_solver(MESH* m)
	{
		Parameters& p = parameters_[m];

		if (p.initialized_ && !p.solver_ready_ && p.selected_free_vertices_set_ &&
			p.selected_free_vertices_set_->size() > 0 && p.selected_handle_vertices_set_ &&
			p.selected_handle_vertices_set_->size() > 0)
		{
			CellMarkerStore<MESH, Vertex> working_vertices_marker(*m);

			// check that handle vertices are surrounded only by handle or free vertices
			bool handle_ok = true;
			foreach_cell(*m, [&](Vertex v) -> bool {
				if (p.selected_handle_vertices_set_->contains(v))
				{
					foreach_adjacent_vertex_through_edge(*m, v, [&](Vertex av) -> bool {
						if (!p.selected_handle_vertices_set_->contains(av) &&
							!p.selected_free_vertices_set_->contains(av))
							handle_ok = false;
						return handle_ok;
					});
				}
				return handle_ok;
			});
			if (!handle_ok)
			{
				std::cout << "surface_deformation: handle is not defined in the free area";
				return;
			}

			// build the cell cache of working area vertices (and mark them)
			p.working_cells_->template build<Vertex>([&](Vertex v) -> bool {
				if (p.selected_handle_vertices_set_->contains(v)) // handle vertices
				{
					working_vertices_marker.mark(v);
					return true;
				}
				if (p.selected_free_vertices_set_->contains(v)) // free vertices
				{
					working_vertices_marker.mark(v);
					foreach_adjacent_vertex_through_edge(
						*m, v,
						[&](Vertex av) -> bool // and their 2-ring
						{
							if (!p.selected_free_vertices_set_->contains(av) &&
								!p.selected_handle_vertices_set_->contains(av) &&
								!working_vertices_marker.is_marked(av))
							{
								p.working_cells_->add(av);
								working_vertices_marker.mark(av);
								foreach_adjacent_vertex_through_edge(*m, av, [&](Vertex aav) -> bool {
									if (!p.selected_free_vertices_set_->contains(aav) &&
										!p.selected_handle_vertices_set_->contains(aav) &&
										!working_vertices_marker.is_marked(aav))
									{
										p.working_cells_->add(aav);
										working_vertices_marker.mark(aav);
									}
									return true;
								});
							}
							return true;
						});
					return true;
				}
				return false;
			});

			// build the cell cache of working area edges
			p.working_cells_->template build<Edge>([&](Edge e) -> bool {
				auto vertices = incident_vertices(*m, e);
				return (working_vertices_marker.is_marked(vertices[0]) &&
						working_vertices_marker.is_marked(vertices[1]));
			});

			// index the working area vertices
			uint32 nb_vertices = 0;
			// start with the free vertices
			foreach_cell(*p.working_cells_, [&](Vertex v) -> bool {
				if (p.selected_free_vertices_set_->contains(v))
					value<uint32>(*m, p.vertex_index_, v) = nb_vertices++;
				return true;
			});
			// then the others (handle & area boundary <=> constrained)
			foreach_cell(*p.working_cells_, [&](Vertex v) -> bool {
				if (!p.selected_free_vertices_set_->contains(v))
					value<uint32>(*m, p.vertex_index_, v) = nb_vertices++;
				return true;
			});

			// init laplacian matrix
			p.working_LAPL_.setZero();
			p.working_LAPL_.resize(nb_vertices, nb_vertices);
			std::vector<Eigen::Triplet<Scalar>> LAPLcoeffs;
			LAPLcoeffs.reserve(nb_vertices * 10);
			foreach_cell(*p.working_cells_, [&](Edge e) -> bool {
				Scalar w = value<Scalar>(*m, p.edge_weight_, e);
				auto vertices = incident_vertices(*m, e);
				uint32 vidx1 = value<uint32>(*m, p.vertex_index_, vertices[0]);
				uint32 vidx2 = value<uint32>(*m, p.vertex_index_, vertices[1]);
				LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx2), w));
				LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx1), w));
				LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx1), int(vidx1), -w));
				LAPLcoeffs.push_back(Eigen::Triplet<Scalar>(int(vidx2), int(vidx2), -w));
				return true;
			});
			p.working_LAPL_.setFromTriplets(LAPLcoeffs.begin(), LAPLcoeffs.end());

			// init bi-laplacian matrix
			p.working_BILAPL_.setZero();
			p.working_BILAPL_.resize(nb_vertices, nb_vertices);
			p.working_BILAPL_ = p.working_LAPL_ * p.working_LAPL_;

			// set constrained vertices
			foreach_cell(*p.working_cells_, [&](Vertex v) -> bool {
				if (!p.selected_free_vertices_set_->contains(v))
				{
					int idx = int(value<uint32>(*m, p.vertex_index_, v));
					p.working_LAPL_.prune([&](int i, int, Scalar) { return i != idx; });
					p.working_LAPL_.coeffRef(idx, idx) = 1.0;
					p.working_BILAPL_.prune([&](int i, int, Scalar) { return i != idx; });
					p.working_BILAPL_.coeffRef(idx, idx) = 1.0;
				}
				return true;
			});

			p.working_LAPL_.makeCompressed();
			p.working_BILAPL_.makeCompressed();

			if (p.solver_)
				delete p.solver_;
			p.solver_ = new Eigen::SparseLU<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>>(p.working_BILAPL_);

			p.solver_ready_ = true;
		}
	}

	void as_rigid_as_possible(MESH* m)
	{
		Parameters& p = parameters_[m];

		if (!p.initialized_)
			return;

		if (!p.solver_ready_)
			build_solver(m);

		if (!p.solver_ready_)
			return;

		parallel_foreach_cell(*p.working_cells_, [&](Vertex v) -> bool {
			Mat3 cov;
			cov.setZero();
			const Vec3& pos = value<Vec3>(*m, p.vertex_position_, v);
			const Vec3& pos_i = value<Vec3>(*m, p.vertex_position_init_, v);
			foreach_adjacent_vertex_through_edge(*m, v, [&](Vertex av) -> bool {
				Vec3 vec = value<Vec3>(*m, p.vertex_position_, av) - pos;
				Vec3 vec_i = value<Vec3>(*m, p.vertex_position_init_, av) - pos_i;
				for (uint32 i = 0; i < 3; ++i)
					for (uint32 j = 0; j < 3; ++j)
						cov(i, j) += vec[i] * vec_i[j];
				return true;
			});
			Eigen::JacobiSVD<Mat3> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
			Mat3 R = svd.matrixU() * svd.matrixV().transpose();
			if (R.determinant() < 0)
			{
				Mat3 U = svd.matrixU();
				for (uint32 i = 0; i < 3; ++i)
					U(i, 2) *= -1;
				R = U * svd.matrixV().transpose();
			}
			value<Mat3>(*m, p.vertex_rotation_matrix_, v) = R;
			return true;
		});

		parallel_foreach_cell(*p.working_cells_, [&](Vertex v) -> bool {
			uint32 degree = 0;
			Mat3 r;
			r.setZero();
			foreach_adjacent_vertex_through_edge(*m, v, [&](Vertex av) -> bool {
				r += value<Mat3>(*m, p.vertex_rotation_matrix_, av);
				++degree;
				return true;
			});
			r += value<Mat3>(*m, p.vertex_rotation_matrix_, v);
			r /= degree + 1;
			value<Vec3>(*m, p.vertex_rotated_diff_coord_, v) = r * value<Vec3>(*m, p.vertex_diff_coord_, v);
			return true;
		});

		uint32 nb_vertices = p.working_cells_->template size<Vertex>();

		Eigen::MatrixXd rdiff(nb_vertices, 3);
		parallel_foreach_cell(*p.working_cells_, [&](Vertex v) -> bool {
			const Vec3& rdcv = value<Vec3>(*m, p.vertex_rotated_diff_coord_, v);
			uint32 vidx = value<uint32>(*m, p.vertex_index_, v);
			rdiff(vidx, 0) = rdcv[0];
			rdiff(vidx, 1) = rdcv[1];
			rdiff(vidx, 2) = rdcv[2];
			return true;
		});
		Eigen::MatrixXd rbdiff(nb_vertices, 3);
		rbdiff = p.working_LAPL_ * rdiff;
		parallel_foreach_cell(*p.working_cells_, [&](Vertex v) -> bool {
			Vec3& rbdcv = value<Vec3>(*m, p.vertex_rotated_bi_diff_coord_, v);
			uint32 vidx = value<uint32>(*m, p.vertex_index_, v);
			rbdcv[0] = rbdiff(vidx, 0);
			rbdcv[1] = rbdiff(vidx, 1);
			rbdcv[2] = rbdiff(vidx, 2);
			return true;
		});

		Eigen::MatrixXd x(nb_vertices, 3);
		Eigen::MatrixXd b(nb_vertices, 3);

		parallel_foreach_cell(*p.working_cells_, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(*m, p.vertex_index_, v);
			if (p.selected_free_vertices_set_->contains(v))
			{
				const Vec3& rbdc = value<Vec3>(*m, p.vertex_rotated_bi_diff_coord_, v);
				b.coeffRef(vidx, 0) = rbdc[0];
				b.coeffRef(vidx, 1) = rbdc[1];
				b.coeffRef(vidx, 2) = rbdc[2];
			}
			else
			{
				const Vec3& pos = value<Vec3>(*m, p.vertex_position_, v);
				b.coeffRef(vidx, 0) = pos[0];
				b.coeffRef(vidx, 1) = pos[1];
				b.coeffRef(vidx, 2) = pos[2];
			}
			return true;
		});

		x = p.solver_->solve(b);

		parallel_foreach_cell(*p.working_cells_, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(*m, p.vertex_index_, v);
			Vec3& pos = value<Vec3>(*m, p.vertex_position_, v);
			pos[0] = x(vidx, 0);
			pos[1] = x(vidx, 1);
			pos[2] = x(vidx, 2);
			return true;
		});
	}

public:
	void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{
		Parameters& p = parameters_[&m];
		p.vertex_position_ = vertex_position;
	}

	void set_selected_free_vertices_set(const MESH& m, CellsSet<MESH, Vertex>* set)
	{
		Parameters& p = parameters_[&m];
		p.selected_free_vertices_set_ = set;
	}

	void set_selected_handle_vertices_set(const MESH& m, CellsSet<MESH, Vertex>* set)
	{
		Parameters& p = parameters_[&m];
		p.selected_handle_vertices_set_ = set;
	}

protected:
	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
		mesh_provider_->foreach_mesh([this](MESH* m, const std::string&) { init_mesh(m); });
		connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
			mesh_provider_, this, &SurfaceDeformation<MESH>::init_mesh));
	}

	void key_press_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_D)
		{
			if (selected_mesh_)
			{
				Parameters& p = parameters_[selected_mesh_];
				if (p.vertex_position_ && p.selected_handle_vertices_set_ &&
					p.selected_handle_vertices_set_->size() > 0)
				{
					drag_z_ = 0.0;
					p.selected_handle_vertices_set_->foreach_cell([&](Vertex v) {
						const Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_, v);
						rendering::GLVec4d vec(pos[0], pos[1], pos[2], 1.0);
						vec = view->projection_matrix_d() * view->modelview_matrix_d() * vec;
						vec /= vec[3];
						drag_z_ += (1.0 + vec[2]) / 2.0;
					});
					drag_z_ /= p.selected_handle_vertices_set_->size();
					previous_drag_pos_ = view->unproject(view->previous_mouse_x(), view->previous_mouse_y(), drag_z_);

					for (View* v : linked_views_)
						v->lock_scene_bb();

					dragging_ = true;
				}
			}
		}
	}

	void key_release_event(View* view, int32 key_code) override
	{
		unused_parameters(view);
		if (key_code == GLFW_KEY_D)
		{
			if (dragging_)
			{
				dragging_ = false;

				for (View* v : linked_views_)
					v->unlock_scene_bb();
			}
		}
	}

	void mouse_move_event(View* view, int32 x, int32 y) override
	{
		if (dragging_)
		{
			Parameters& p = parameters_[selected_mesh_];

			rendering::GLVec3d drag_pos = view->unproject(x, y, drag_z_);
			Vec3 t = drag_pos - previous_drag_pos_;
			p.selected_handle_vertices_set_->foreach_cell(
				[&](Vertex v) { value<Vec3>(*selected_mesh_, p.vertex_position_, v) += t; });
			as_rigid_as_possible(selected_mesh_);
			previous_drag_pos_ = drag_pos;

			mesh_provider_->emit_attribute_changed(selected_mesh_, p.vertex_position_.get());
		}
	}

	void interface() override
	{
		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});

		if (ImGui::ListBoxHeader("Mesh"))
		{
			mesh_provider_->foreach_mesh([this](MESH* m, const std::string& name) {
				if (ImGui::Selectable(name.c_str(), m == selected_mesh_))
					selected_mesh_ = m;
			});
			ImGui::ListBoxFooter();
		}

		if (selected_mesh_)
		{
			float X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

			MeshData<MESH>* md = mesh_provider_->mesh_data(selected_mesh_);
			Parameters& p = parameters_[selected_mesh_];

			if (ImGui::BeginCombo("Position", p.vertex_position_ ? p.vertex_position_->name().c_str() : "-- select --"))
			{
				foreach_attribute<Vec3, Vertex>(*selected_mesh_,
												[&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
													bool is_selected = attribute == p.vertex_position_;
													if (ImGui::Selectable(attribute->name().c_str(), is_selected))
														set_vertex_position(*selected_mesh_, attribute);
													if (is_selected)
														ImGui::SetItemDefaultFocus();
												});
				ImGui::EndCombo();
			}
			if (p.vertex_position_)
			{
				ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
				if (ImGui::Button("X##position"))
					set_vertex_position(*selected_mesh_, nullptr);

				ImGui::Separator();

				if (ImGui::Button("Initialize"))
					initialize_mesh_data(selected_mesh_);

				ImGui::Separator();

				if (ImGui::BeginCombo("Free vertices", p.selected_free_vertices_set_
														   ? p.selected_free_vertices_set_->name().c_str()
														   : "-- select --"))
				{
					md->template foreach_cells_set<Vertex>([&](CellsSet<MESH, Vertex>& cs) {
						bool is_selected = &cs == p.selected_free_vertices_set_;
						if (ImGui::Selectable(cs.name().c_str(), is_selected))
							set_selected_free_vertices_set(*selected_mesh_, &cs);
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					});
					ImGui::EndCombo();
				}
				if (p.selected_free_vertices_set_)
				{
					ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
					if (ImGui::Button("X##selected_free_vertices_set"))
						set_selected_free_vertices_set(*selected_mesh_, nullptr);
				}

				if (ImGui::BeginCombo("Handle vertices", p.selected_handle_vertices_set_
															 ? p.selected_handle_vertices_set_->name().c_str()
															 : "-- select --"))
				{
					md->template foreach_cells_set<Vertex>([&](CellsSet<MESH, Vertex>& cs) {
						bool is_selected = &cs == p.selected_handle_vertices_set_;
						if (ImGui::Selectable(cs.name().c_str(), is_selected))
							set_selected_handle_vertices_set(*selected_mesh_, &cs);
						if (is_selected)
							ImGui::SetItemDefaultFocus();
					});
					ImGui::EndCombo();
				}
				if (p.selected_handle_vertices_set_)
				{
					ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
					if (ImGui::Button("X##selected_handle_vertices_set"))
						set_selected_handle_vertices_set(*selected_mesh_, nullptr);
					ImGui::TextUnformatted("Press D to drag the handle");
				}
			}
		}

		ImGui::End();
	}

private:
	MESH* selected_mesh_;
	std::unordered_map<const MESH*, Parameters> parameters_;
	std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
	MeshProvider<MESH>* mesh_provider_;

	bool dragging_;
	float64 drag_z_;
	rendering::GLVec3d previous_drag_pos_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SURFACE_DEFORMATION_H_

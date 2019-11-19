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

#ifndef CGOGN_MODULE_SHALLOW_WATER_H_
#define CGOGN_MODULE_SHALLOW_WATER_H_

#include <cgogn/ui/module.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>
#include <cgogn/core/functions/traversals/halfedge.h>
#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/length.h>

// #include <cgogn/core/utils/timer.h>

#include <cgogn/simulation/algos/shallow_water/riemann_solver.h>

namespace cgogn
{

namespace ui
{

template <typename MESH>
class ShallowWater : public ViewModule
{
	static_assert(mesh_traits<MESH>::dimension == 2, "ShallowWater can only be used with meshes of dimension 2");

    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

    using Vertex = typename mesh_traits<MESH>::Vertex;
    using HalfEdge = typename mesh_traits<MESH>::HalfEdge;
    using Edge = typename mesh_traits<MESH>::Edge;
    using Face = typename mesh_traits<MESH>::Face;

    using Vec3 = geometry::Vec3;
    using Scalar = geometry::Scalar;

	using BoundaryCondition = simulation::BoundaryCondition;

public:

	ShallowWater(const App& app) :
		ViewModule(app, "ShallowWater (" + std::string{mesh_traits<MESH>::name} + ")")
	{}
	~ShallowWater()
	{}

protected:

	void init() override
	{
		mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));

		domain2D_ = mesh_provider_->load_surface_from_file("/Users/kraemer/Data/surface/grid_tri_high.off");

		vertex_position_ = get_attribute<Vec3, Vertex>(*domain2D_, "position");
		vertex_water_position_ = add_attribute<Vec3, Vertex>(*domain2D_, "water_position");

		face_phi_ = add_attribute<Scalar, Face>(*domain2D_, "phi");
		face_zb_ = add_attribute<Scalar, Face>(*domain2D_, "zb");
		face_h_ = add_attribute<Scalar, Face>(*domain2D_, "h");
		face_q_ = add_attribute<Scalar, Face>(*domain2D_, "q");
		face_r_ = add_attribute<Scalar, Face>(*domain2D_, "r");
		face_centroid_ = add_attribute<Vec3, Face>(*domain2D_, "centroid");
		face_area_ = add_attribute<Scalar, Face>(*domain2D_, "area");
		face_swept_ = add_attribute<Scalar, Face>(*domain2D_, "swept");
		face_discharge_ = add_attribute<Scalar, Face>(*domain2D_, "discharge");

		edge_f1_ = add_attribute<Scalar, Edge>(*domain2D_, "f1");
		edge_f2_ = add_attribute<Scalar, Edge>(*domain2D_, "f2");
		edge_f3_ = add_attribute<Scalar, Edge>(*domain2D_, "f3");
		edge_s2L_ = add_attribute<Scalar, Edge>(*domain2D_, "s2L");
		edge_s2R_ = add_attribute<Scalar, Edge>(*domain2D_, "s2R");
		edge_normX_ = add_attribute<Scalar, Edge>(*domain2D_, "normX");
		edge_normY_ = add_attribute<Scalar, Edge>(*domain2D_, "normY");
		edge_length_ = add_attribute<Scalar, Edge>(*domain2D_, "length");
		edge_bc_value_ = add_attribute<Scalar, Edge>(*domain2D_, "bc_value");
		edge_bc_type_ = add_attribute<BoundaryCondition, Edge>(*domain2D_, "bc_type");
		edge_left_face_ = add_attribute<uint32, Edge>(*domain2D_, "left_face");

		vertex_water_position_->copy(vertex_position_.get());

		parallel_foreach_cell(*domain2D_, [&] (Edge e) -> bool
		{
			if (is_incident_to_boundary(*domain2D_, e))
				value<BoundaryCondition>(*domain2D_, edge_bc_type_, e) = simulation::BC_Q;
			
			Scalar l = geometry::length(*domain2D_, e, vertex_position_.get());
			std::vector<Vertex> vertices = incident_vertices(*domain2D_, e);
			Vec3 vec = value<Vec3>(*domain2D_, vertex_position_, vertices[1]) - value<Vec3>(*domain2D_, vertex_position_, vertices[0]);
			value<Scalar>(*domain2D_, edge_length_, e) = l;
			value<Scalar>(*domain2D_, edge_normX_, e) = vec[1] / l;
			value<Scalar>(*domain2D_, edge_normY_, e) = -vec[0] / l;
			value<uint32>(*domain2D_, edge_left_face_, e) = index_of(*domain2D_, incident_faces(*domain2D_, e)[0]);

			return true;
		});

		srand(time(0));
		parallel_foreach_cell(*domain2D_, [&] (Face f) -> bool
		{
			uint32 nbv = 0;
			Scalar zb = 0.0;
			Vec3 centroid{0, 0, 0};
			foreach_incident_vertex(*domain2D_, f, [&] (Vertex v) -> bool
			{
				const Vec3& p = value<Vec3>(*domain2D_, vertex_position_, v);
				zb += p[2];
				centroid += p;
				nbv++;
				return true;
			});
			zb /= nbv;
			centroid /= nbv;

			Scalar h;
			bool border = false;
			foreach_incident_vertex(*domain2D_, f, [&] (Vertex v) -> bool
			{
				if (is_incident_to_boundary(*domain2D_, v))
					border = true;
				return !border;
			});
			if (border) h = static_cast<Scalar>(rand()) / static_cast<Scalar>(RAND_MAX);
			else h = 0.05;
			
			value<Scalar>(*domain2D_, face_phi_, f) = phi_default_;
			value<Scalar>(*domain2D_, face_zb_, f) = zb;
			value<Scalar>(*domain2D_, face_h_, f) = h;
			value<Scalar>(*domain2D_, face_q_, f) = 0.0;
			value<Scalar>(*domain2D_, face_r_, f) = 0.0;
			value<Vec3>(*domain2D_, face_centroid_, f) = centroid;
			value<Scalar>(*domain2D_, face_area_, f) = geometry::area(*domain2D_, f, vertex_position_.get());

			return true;
		});

		update_render_data();

		running_ = false;

		timer_connection_ = boost::synapse::connect<App::timer_tick>(&app_, [this] ()
		{
			update_render_data();
			for (View* v : linked_views_)
				v->request_update();
		});
	}

	void start()
	{
		running_ = true;

		launch_thread([this] ()
		{
			while (this->running_)
				this->execute_time_step();
		});

		app_.start_timer(50, [this] () -> bool { return !running_; });
	}

	void stop()
	{
		running_ = false;
	}

	void update_time_step()
	{
		face_discharge_->fill(0.);
		face_swept_->fill(0.);

		parallel_foreach_cell(*domain2D_, [&] (Face f) -> bool
		{
			uint32 fidx = index_of(*domain2D_, f);
			foreach_incident_edge(*domain2D_, f, [&] (Edge ie) -> bool
			{
				uint32 ieidx = index_of(*domain2D_, ie);
				Scalar le = (*edge_length_)[ieidx];
				Scalar f1e = (*edge_f1_)[ieidx];
				Scalar lambda = 0.;
				Scalar h = (*face_h_)[fidx];
				if (h > hmin_)
					lambda = fabs((*face_q_)[fidx]*(*edge_normX_)[ieidx] + (*face_r_)[fidx]*(*edge_normY_)[ieidx]) / std::max(h, hmin_) + sqrt(9.81*h);
				(*face_swept_)[fidx] += le * lambda;
				if (fidx == (*edge_left_face_)[ieidx])
					(*face_discharge_)[fidx] -= le * f1e;
				else
					(*face_discharge_)[fidx] += le * f1e;
				return true;
			});
			return true;
		});

		std::vector<Scalar> min_dt_per_thread(thread_pool()->nb_workers());
		for (Scalar& d : min_dt_per_thread) d = std::min(dt_max_, t_max_ - t_); // Timestep for ending simulation

		parallel_foreach_cell(*domain2D_, [&] (Face f) -> bool
		{
			uint32 fidx = index_of(*domain2D_, f);
			uint32 threadidx = current_worker_index();
			// Ensure CFL condition
			Scalar cfl = (*face_area_)[fidx] / std::max((*face_swept_)[fidx], small_);
			min_dt_per_thread[threadidx] = std::min(min_dt_per_thread[threadidx], cfl);
			// Ensure overdry condition
			if ((*face_area_)[fidx]*(*face_phi_)[fidx]*((*face_h_)[fidx]+(*face_zb_)[fidx]) < (-(*face_discharge_)[fidx]*min_dt_per_thread[threadidx]))
				min_dt_per_thread[threadidx] = - (*face_area_)[fidx]*(*face_phi_)[fidx]*((*face_h_)[fidx]+(*face_zb_)[fidx]) / (*face_discharge_)[fidx];
			return true;
		});

		dt_ = *(std::min_element(min_dt_per_thread.begin(), min_dt_per_thread.end()));
	}

	void execute_time_step()
	{
		auto start = std::chrono::high_resolution_clock::now();

		parallel_foreach_cell(*domain2D_, [&] (Edge e) -> bool
		{
			uint32 eidx = index_of(*domain2D_, e);

			// solve flux on edge
			simulation::Str_Riemann_Flux riemann_flux;

			if (is_incident_to_boundary(*domain2D_, e)) // border conditions
			{
				uint32 fidx = (*edge_left_face_)[eidx];
				if ((*face_phi_)[fidx] > small_)
					riemann_flux = simulation::border_condition((*edge_bc_type_)[eidx], (*edge_bc_value_)[eidx], (*edge_normX_)[eidx], (*edge_normY_)[eidx], (*face_q_)[fidx], (*face_r_)[fidx], (*face_h_)[fidx]+(*face_zb_)[fidx], (*face_zb_)[fidx], 9.81, hmin_, small_);
			}
			else // Inner cell: use the lateralised Riemann solver
			{
				std::vector<Face> faces = incident_faces(*domain2D_, e);
				uint32 f1idx = index_of(*domain2D_, faces[0]);
				uint32 f2idx = index_of(*domain2D_, faces[1]);

				Scalar phiL = (*face_phi_)[f1idx];
				Scalar phiR = (*face_phi_)[f2idx];
				Scalar zbL = (*face_zb_)[f1idx];
				Scalar zbR = (*face_zb_)[f2idx];
				if ((*face_h_)[f1idx] > hmin_ || (*face_h_)[f2idx] > hmin_)
				{
					Scalar hL = (*face_h_)[f1idx];
					Scalar hR = (*face_h_)[f2idx];
					Scalar qL = (*face_q_)[f1idx]*(*edge_normX_)[eidx] + (*face_r_)[f1idx]*(*edge_normY_)[eidx];
					Scalar qR = (*face_q_)[f2idx]*(*edge_normX_)[eidx] + (*face_r_)[f2idx]*(*edge_normY_)[eidx];
					Scalar rL = -(*face_q_)[f1idx]*(*edge_normY_)[eidx] + (*face_r_)[f1idx]*(*edge_normX_)[eidx];
					Scalar rR = -(*face_q_)[f2idx]*(*edge_normY_)[eidx] + (*face_r_)[f2idx]*(*edge_normX_)[eidx];

					riemann_flux = simulation::Solv_HLLC(9.81, hmin_, small_, zbL, zbR, phiL, phiR, hL, qL, rL, hR, qR, rR);
				}
			}

			(*edge_f1_)[eidx] = riemann_flux.F1;
			(*edge_f2_)[eidx] = riemann_flux.F2;
			(*edge_f3_)[eidx] = riemann_flux.F3;
			(*edge_s2L_)[eidx] = riemann_flux.s2L;
			(*edge_s2R_)[eidx] = riemann_flux.s2R;

			return true;
		});

		update_time_step();

		// simu_data_access_.lock();

		parallel_foreach_cell(*domain2D_, [&] (Face f) -> bool
		{
			uint32 fidx = index_of(*domain2D_, f);

			foreach_incident_edge(*domain2D_, f, [&] (Edge ie) -> bool
			{
				uint32 ieidx = index_of(*domain2D_, ie);
				Scalar fact = dt_ * (*edge_length_)[ieidx];
				Scalar factF = 0.;
				if ((*face_phi_)[fidx] > small_)
					factF = fact / (*face_area_)[fidx] * (*face_phi_)[fidx];
				if (fidx == (*edge_left_face_)[ieidx])
				{
					(*face_h_)[fidx] -= factF * (*edge_f1_)[ieidx];
					(*face_q_)[fidx] -= factF * (((*edge_f2_)[ieidx] + (*edge_s2L_)[ieidx])*(*edge_normX_)[ieidx] - (*edge_f3_)[ieidx]*(*edge_normY_)[ieidx]);
					(*face_r_)[fidx] -= factF * ((*edge_f3_)[ieidx]*(*edge_normX_)[ieidx] + ((*edge_f2_)[ieidx] + (*edge_s2L_)[ieidx])*(*edge_normY_)[ieidx]);
				}
				else
				{
					(*face_h_)[fidx] += factF * (*edge_f1_)[ieidx];
					(*face_q_)[fidx] += factF * (( (*edge_f2_)[ieidx] + (*edge_s2R_)[ieidx])*(*edge_normX_)[ieidx] - (*edge_f3_)[ieidx]*(*edge_normY_)[ieidx]);
					(*face_r_)[fidx] += factF * ( (*edge_f3_)[ieidx]*(*edge_normX_)[ieidx] + ( (*edge_f2_)[ieidx]+(*edge_s2R_)[ieidx])*(*edge_normY_)[ieidx]);
				}
				return true;
			});
			return true;
		});

		parallel_foreach_cell(*domain2D_, [&] (Face f) -> bool
		{
			uint32 fidx = index_of(*domain2D_, f);

			// friction
			if (friction_ != 0)
			{
				Scalar qx = (*face_q_)[fidx]*cos(alphaK_) + (*face_r_)[fidx]*sin(alphaK_);
				Scalar qy = - (*face_q_)[fidx]*sin(alphaK_) + (*face_r_)[fidx]*cos(alphaK_);
				if ((*face_h_)[fidx] > hmin_)
				{
					qx = qx * exp(-(9.81 * sqrt(qx*qx+qy*qy) / (std::max(kx_*kx_,small_*small_) * pow((*face_h_)[fidx],7./3.))) * dt_);
					qy = qy * exp(-(9.81 * sqrt(qx*qx+qy*qy) / (std::max(ky_*ky_,small_*small_) * pow((*face_h_)[fidx],7./3.))) * dt_);
				}
				else
				{
					qx = 0.;
					qy = 0.;
				}
				(*face_q_)[fidx] = qx*cos(alphaK_) - qy*sin(alphaK_);
				(*face_r_)[fidx] = qx*sin(alphaK_) + qy*cos(alphaK_);
			}

			// optional correction
			// Negative water depth
			if ((*face_h_)[fidx] < 0.)
			{
				(*face_h_)[fidx] = 0.;
				(*face_q_)[fidx] = 0.;
				(*face_r_)[fidx] = 0.;
			}

			// Abnormal large velocity => Correction of q and r to respect Vmax and Frmax
			if ((*face_h_)[fidx] > hmin_)
			{
				Scalar v = sqrt((*face_q_)[fidx]*(*face_q_)[fidx]+(*face_r_)[fidx]*(*face_r_)[fidx]) / std::max((*face_h_)[fidx], small_);
				Scalar c = sqrt(9.81 * std::max((*face_h_)[fidx], small_));
				Scalar Fr = v / c;
				Scalar Fact = std::max({ 1e0, v / v_max_, Fr / Fr_max_ });
				(*face_q_)[fidx] /= Fact;
				(*face_r_)[fidx] /= Fact;
			}
			else // Quasi-zero
			{
				(*face_q_)[fidx] = 0.;
				(*face_r_)[fidx] = 0.;
			}
			return true;
		});

		// simu_data_access_.unlock();

		// if (adaptive_mesh_)
		// {
		// 	if (nb_iter_ % iteradapt_ == 0 )
		// 	{
		// 		map_->lock_topo_access();
		// 		try_simplification();
		// 		try_subdivision();
		// 		map_->unlock_topo_access();
		// 	}
		// }

		t_ += dt_;
		// nb_iter_++;

		auto end = std::chrono::high_resolution_clock::now();

		std::chrono::nanoseconds sleep_duration =
			std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::duration<Scalar>(dt_))
			- std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

		if (sleep_duration > std::chrono::nanoseconds::zero())
			std::this_thread::sleep_for(sleep_duration);
		
		if (t_ == t_max_)
			stop();
	}

	void update_render_data()
	{
		foreach_cell(*domain2D_, [&] (Vertex v) -> bool
		{
			Scalar h = 0.0;
			uint32 nbf = 0;
			foreach_incident_face(*domain2D_, v, [&] (Face f) -> bool
			{
				h += value<Scalar>(*domain2D_, face_h_, f);
				++nbf;
				return true;
			});
			value<Vec3>(*domain2D_, vertex_water_position_, v)[2] = h / nbf;
			return true;
		});
		
		mesh_provider_->emit_attribute_changed(domain2D_, vertex_water_position_.get());
	}

    void interface() override
	{
		ImGui::Begin(name_.c_str(), nullptr, ImGuiWindowFlags_NoSavedSettings);
		ImGui::SetWindowSize({0, 0});

		if (!running_)
		{
			if (ImGui::Button("Start"))
				start();
			if (ImGui::Button("step"))
			{
				execute_time_step();
				update_render_data();
				for (View* v : linked_views_)
					v->request_update();
			}
		}
		else
		{
			if (ImGui::Button("Stop"))
				stop();
		}
		
		ImGui::End();
	}

private:

	MeshProvider<MESH>* mesh_provider_;
	MESH* domain2D_;

	std::shared_ptr<Attribute<Vec3>> vertex_position_;
	std::shared_ptr<Attribute<Vec3>> vertex_water_position_;

	std::shared_ptr<Attribute<Scalar>> face_phi_;
	std::shared_ptr<Attribute<Scalar>> face_zb_;
	std::shared_ptr<Attribute<Scalar>> face_h_;
	std::shared_ptr<Attribute<Scalar>> face_q_;
	std::shared_ptr<Attribute<Scalar>> face_r_;
	std::shared_ptr<Attribute<Vec3>> face_centroid_;
	std::shared_ptr<Attribute<Scalar>> face_area_;
	std::shared_ptr<Attribute<Scalar>> face_swept_;
	std::shared_ptr<Attribute<Scalar>> face_discharge_;

	std::shared_ptr<Attribute<Scalar>> edge_f1_;
	std::shared_ptr<Attribute<Scalar>> edge_f2_;
	std::shared_ptr<Attribute<Scalar>> edge_f3_;
	std::shared_ptr<Attribute<Scalar>> edge_s2L_;
	std::shared_ptr<Attribute<Scalar>> edge_s2R_;
	std::shared_ptr<Attribute<Scalar>> edge_normX_;
	std::shared_ptr<Attribute<Scalar>> edge_normY_;
	std::shared_ptr<Attribute<Scalar>> edge_length_;
	std::shared_ptr<Attribute<Scalar>> edge_bc_value_;
	std::shared_ptr<Attribute<BoundaryCondition>> edge_bc_type_;
	std::shared_ptr<Attribute<uint32>> edge_left_face_;

	Scalar phi_default_ = 1.0;
	Scalar kx_ = 35.0;
	Scalar ky_ = 35.0;
	Scalar alphaK_ = 0.0;
	Scalar hmin_ = 1e-3;
	Scalar small_ = 1e-35;
	uint32 friction_ = 0;
	Scalar v_max_ = 10.0;
	Scalar Fr_max_ = 5.0;

	Scalar t_ = 0.0;
	Scalar dt_;
	Scalar t_max_ = 5.0;
	Scalar dt_max_ = 1.0;

	std::shared_ptr<boost::synapse::connection> timer_connection_;
	bool running_ = false;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SHALLOW_WATER_H_

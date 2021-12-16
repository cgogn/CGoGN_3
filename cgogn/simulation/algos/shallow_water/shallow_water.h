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

#ifndef CGOGN_SIMULATION_SHALLOW_WATER_H_
#define CGOGN_SIMULATION_SHALLOW_WATER_H_

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/mesh_info.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/length.h>

#include <cgogn/simulation/algos/shallow_water/riemann_solver.h>

namespace cgogn
{

namespace simulation
{

namespace shallow_water
{

using Vec3 = geometry::Vec3;
using Scalar = geometry::Scalar;

template <typename MESH>
struct Attributes
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;

	std::shared_ptr<Attribute<Vec3>> vertex_position_;

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
	std::shared_ptr<Attribute<uint32>> edge_left_face_index_;
};

struct Context
{
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
	Scalar t_max_ = 10.0;
	Scalar dt_ = 0.0;
	Scalar dt_max_ = 1.0;
};

template <typename MESH>
void get_attributes(MESH& m, Attributes<MESH>& swa)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	swa.vertex_position_ = get_attribute<Vec3, Vertex>(m, "position");

	swa.face_phi_ = add_attribute<Scalar, Face>(m, "phi");
	swa.face_zb_ = add_attribute<Scalar, Face>(m, "zb");
	swa.face_h_ = add_attribute<Scalar, Face>(m, "h");
	swa.face_q_ = add_attribute<Scalar, Face>(m, "q");
	swa.face_r_ = add_attribute<Scalar, Face>(m, "r");
	swa.face_centroid_ = add_attribute<Vec3, Face>(m, "centroid");
	swa.face_area_ = add_attribute<Scalar, Face>(m, "area");
	swa.face_swept_ = add_attribute<Scalar, Face>(m, "swept");
	swa.face_discharge_ = add_attribute<Scalar, Face>(m, "discharge");

	swa.edge_f1_ = add_attribute<Scalar, Edge>(m, "f1");
	swa.edge_f2_ = add_attribute<Scalar, Edge>(m, "f2");
	swa.edge_f3_ = add_attribute<Scalar, Edge>(m, "f3");
	swa.edge_s2L_ = add_attribute<Scalar, Edge>(m, "s2L");
	swa.edge_s2R_ = add_attribute<Scalar, Edge>(m, "s2R");
	swa.edge_normX_ = add_attribute<Scalar, Edge>(m, "normX");
	swa.edge_normY_ = add_attribute<Scalar, Edge>(m, "normY");
	swa.edge_length_ = add_attribute<Scalar, Edge>(m, "length");
	swa.edge_bc_value_ = add_attribute<Scalar, Edge>(m, "bc_value");
	swa.edge_bc_type_ = add_attribute<BoundaryCondition, Edge>(m, "bc_type");
	swa.edge_left_face_index_ = add_attribute<uint32, Edge>(m, "left_face_index");
}

template <typename MESH>
void init_attributes(MESH& m, Attributes<MESH>& swa, Context& swc)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	parallel_foreach_cell(m, [&](Edge e) -> bool {
		if (is_incident_to_boundary(m, e))
			value<BoundaryCondition>(m, swa.edge_bc_type_, e) = BC_Q;

		std::vector<Vertex> vertices = incident_vertices(m, e);
		Vec3 vec =
			value<Vec3>(m, swa.vertex_position_, vertices[1]) - value<Vec3>(m, swa.vertex_position_, vertices[0]);
		Scalar l = vec.norm();
		value<Scalar>(m, swa.edge_length_, e) = l;
		value<Scalar>(m, swa.edge_normX_, e) = vec[1] / l;
		value<Scalar>(m, swa.edge_normY_, e) = -vec[0] / l;
		value<uint32>(m, swa.edge_left_face_index_, e) = index_of(m, incident_faces(m, e)[0]);

		return true;
	});

	srand(time(0));
	parallel_foreach_cell(m, [&](Face f) -> bool {
		value<Scalar>(m, swa.face_phi_, f) = swc.phi_default_;

		uint32 nbv = 0;
		Scalar zb = 0.0;
		Vec3 centroid{0, 0, 0};
		foreach_incident_vertex(m, f, [&](Vertex v) -> bool {
			const Vec3& p = value<Vec3>(m, swa.vertex_position_, v);
			zb += p[2];
			centroid += p;
			nbv++;
			return true;
		});
		zb /= nbv;
		centroid /= nbv;

		value<Scalar>(m, swa.face_zb_, f) = zb;
		value<Vec3>(m, swa.face_centroid_, f) = centroid;
		value<Scalar>(m, swa.face_area_, f) = geometry::area(m, f, swa.vertex_position_.get());

		Scalar h;
		bool border = false;
		foreach_incident_vertex(m, f, [&](Vertex v) -> bool {
			if (is_incident_to_boundary(m, v))
				border = true;
			return !border;
		});
		if (border)
			h = static_cast<Scalar>(rand()) / static_cast<Scalar>(RAND_MAX) * 2.0;
		else
			h = 0.25;

		value<Scalar>(m, swa.face_h_, f) = h;
		value<Scalar>(m, swa.face_q_, f) = 0.0;
		value<Scalar>(m, swa.face_r_, f) = 0.0;

		return true;
	});
}

template <typename MESH>
void domain_geometry_changed(MESH& m, Attributes<MESH>& swa, Context& swc)
{
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	parallel_foreach_cell(m, [&](Edge e) -> bool {
		std::vector<Vertex> vertices = incident_vertices(m, e);
		Vec3 vec =
			value<Vec3>(m, swa.vertex_position_, vertices[1]) - value<Vec3>(m, swa.vertex_position_, vertices[0]);
		Scalar l = vec.norm();
		value<Scalar>(m, swa.edge_length_, e) = l;
		value<Scalar>(m, swa.edge_normX_, e) = vec[1] / l;
		value<Scalar>(m, swa.edge_normY_, e) = -vec[0] / l;
		return true;
	});

	parallel_foreach_cell(m, [&](Face f) -> bool {
		uint32 nbv = 0;
		Scalar zb = 0.0;
		Vec3 centroid{0, 0, 0};
		foreach_incident_vertex(m, f, [&](Vertex v) -> bool {
			const Vec3& p = value<Vec3>(m, swa.vertex_position_, v);
			zb += p[2];
			centroid += p;
			nbv++;
			return true;
		});
		zb /= nbv;
		centroid /= nbv;

		value<Scalar>(m, swa.face_zb_, f) = zb;
		value<Vec3>(m, swa.face_centroid_, f) = centroid;

		Scalar old_area = value<Scalar>(m, swa.face_area_, f);
		Scalar old_h = value<Scalar>(m, swa.face_h_, f);

		Scalar new_area = std::max(geometry::area(m, f, swa.vertex_position_.get()), swc.small_);
		Scalar new_h = std::max(old_area / new_area * old_h, swc.hmin_);

		value<Scalar>(m, swa.face_area_, f) = new_area;
		value<Scalar>(m, swa.face_h_, f) = new_h;

		return true;
	});
}

template <typename MESH>
void update_time_step(MESH& m, Attributes<MESH>& swa, Context& swc)
{
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	swa.face_discharge_->fill(0.);
	swa.face_swept_->fill(0.);

	parallel_foreach_cell(m, [&](Face f) -> bool {
		uint32 fidx = index_of(m, f);
		foreach_incident_edge(m, f, [&](Edge ie) -> bool {
			uint32 ieidx = index_of(m, ie);
			Scalar le = (*swa.edge_length_)[ieidx];
			Scalar f1e = (*swa.edge_f1_)[ieidx];
			Scalar lambda = 0.;
			Scalar h = (*swa.face_h_)[fidx];
			if (h > swc.hmin_)
				lambda = fabs((*swa.face_q_)[fidx] * (*swa.edge_normX_)[ieidx] +
							  (*swa.face_r_)[fidx] * (*swa.edge_normY_)[ieidx]) /
							 std::max(h, swc.hmin_) +
						 sqrt(9.81 * h);
			(*swa.face_swept_)[fidx] += le * lambda;
			if (fidx == (*swa.edge_left_face_index_)[ieidx])
				(*swa.face_discharge_)[fidx] -= le * f1e;
			else
				(*swa.face_discharge_)[fidx] += le * f1e;
			return true;
		});
		return true;
	});

	std::vector<Scalar> min_dt_per_thread(thread_pool()->nb_workers());
	// for (Scalar& d : min_dt_per_thread) d = std::min(swc.dt_max_, swc.t_max_ - swc.t_); // Timestep for ending
	// simulation
	for (Scalar& d : min_dt_per_thread)
		d = swc.dt_max_;

	parallel_foreach_cell(m, [&](Face f) -> bool {
		uint32 fidx = index_of(m, f);
		uint32 threadidx = current_worker_index();
		// Ensure CFL condition
		Scalar cfl = (*swa.face_area_)[fidx] / std::max((*swa.face_swept_)[fidx], swc.small_);
		min_dt_per_thread[threadidx] = std::min(min_dt_per_thread[threadidx], cfl);
		// Ensure overdry condition
		if ((*swa.face_area_)[fidx] * (*swa.face_phi_)[fidx] * ((*swa.face_h_)[fidx] + (*swa.face_zb_)[fidx]) <
			(-(*swa.face_discharge_)[fidx] * min_dt_per_thread[threadidx]))
			min_dt_per_thread[threadidx] = -(*swa.face_area_)[fidx] * (*swa.face_phi_)[fidx] *
										   ((*swa.face_h_)[fidx] + (*swa.face_zb_)[fidx]) /
										   (*swa.face_discharge_)[fidx];
		return true;
	});

	swc.dt_ = *(std::min_element(min_dt_per_thread.begin(), min_dt_per_thread.end()));
}

template <typename MESH>
void execute_time_step(MESH& m, Attributes<MESH>& swa, Context& swc)
{
	using Edge = typename mesh_traits<MESH>::Edge;
	using Face = typename mesh_traits<MESH>::Face;

	auto start = std::chrono::high_resolution_clock::now();

	parallel_foreach_cell(m, [&](Edge e) -> bool {
		uint32 eidx = index_of(m, e);

		// solve flux on edge
		Str_Riemann_Flux riemann_flux;

		if (is_incident_to_boundary(m, e)) // border conditions
		{
			uint32 fidx = (*swa.edge_left_face_index_)[eidx];
			if ((*swa.face_phi_)[fidx] > swc.small_)
				riemann_flux = border_condition(
					(*swa.edge_bc_type_)[eidx], (*swa.edge_bc_value_)[eidx], (*swa.edge_normX_)[eidx],
					(*swa.edge_normY_)[eidx], (*swa.face_q_)[fidx], (*swa.face_r_)[fidx],
					(*swa.face_h_)[fidx] + (*swa.face_zb_)[fidx], (*swa.face_zb_)[fidx], 9.81, swc.hmin_, swc.small_);
		}
		else // Inner cell: use the lateralised Riemann solver
		{
			std::vector<Face> faces = incident_faces(m, e);
			uint32 f1idx = index_of(m, faces[0]);
			uint32 f2idx = index_of(m, faces[1]);

			Scalar phiL = (*swa.face_phi_)[f1idx];
			Scalar phiR = (*swa.face_phi_)[f2idx];
			Scalar zbL = (*swa.face_zb_)[f1idx];
			Scalar zbR = (*swa.face_zb_)[f2idx];
			if ((*swa.face_h_)[f1idx] > swc.hmin_ || (*swa.face_h_)[f2idx] > swc.hmin_)
			{
				Scalar hL = (*swa.face_h_)[f1idx];
				Scalar hR = (*swa.face_h_)[f2idx];
				Scalar qL =
					(*swa.face_q_)[f1idx] * (*swa.edge_normX_)[eidx] + (*swa.face_r_)[f1idx] * (*swa.edge_normY_)[eidx];
				Scalar qR =
					(*swa.face_q_)[f2idx] * (*swa.edge_normX_)[eidx] + (*swa.face_r_)[f2idx] * (*swa.edge_normY_)[eidx];
				Scalar rL = -(*swa.face_q_)[f1idx] * (*swa.edge_normY_)[eidx] +
							(*swa.face_r_)[f1idx] * (*swa.edge_normX_)[eidx];
				Scalar rR = -(*swa.face_q_)[f2idx] * (*swa.edge_normY_)[eidx] +
							(*swa.face_r_)[f2idx] * (*swa.edge_normX_)[eidx];

				riemann_flux = Solv_HLLC(9.81, swc.hmin_, swc.small_, zbL, zbR, phiL, phiR, hL, qL, rL, hR, qR, rR);
			}
		}

		(*swa.edge_f1_)[eidx] = riemann_flux.F1;
		(*swa.edge_f2_)[eidx] = riemann_flux.F2;
		(*swa.edge_f3_)[eidx] = riemann_flux.F3;
		(*swa.edge_s2L_)[eidx] = riemann_flux.s2L;
		(*swa.edge_s2R_)[eidx] = riemann_flux.s2R;

		return true;
	});

	update_time_step(m, swa, swc);

	// simu_data_access_.lock();

	parallel_foreach_cell(m, [&](Face f) -> bool {
		uint32 fidx = index_of(m, f);

		foreach_incident_edge(m, f, [&](Edge ie) -> bool {
			uint32 ieidx = index_of(m, ie);
			Scalar fact = swc.dt_ * (*swa.edge_length_)[ieidx];
			Scalar factF = 0.;
			if ((*swa.face_phi_)[fidx] > swc.small_)
				factF = fact / (*swa.face_area_)[fidx] * (*swa.face_phi_)[fidx];
			if (fidx == (*swa.edge_left_face_index_)[ieidx])
			{
				(*swa.face_h_)[fidx] -= factF * (*swa.edge_f1_)[ieidx];
				(*swa.face_q_)[fidx] -=
					factF * (((*swa.edge_f2_)[ieidx] + (*swa.edge_s2L_)[ieidx]) * (*swa.edge_normX_)[ieidx] -
							 (*swa.edge_f3_)[ieidx] * (*swa.edge_normY_)[ieidx]);
				(*swa.face_r_)[fidx] -=
					factF * ((*swa.edge_f3_)[ieidx] * (*swa.edge_normX_)[ieidx] +
							 ((*swa.edge_f2_)[ieidx] + (*swa.edge_s2L_)[ieidx]) * (*swa.edge_normY_)[ieidx]);
			}
			else
			{
				(*swa.face_h_)[fidx] += factF * (*swa.edge_f1_)[ieidx];
				(*swa.face_q_)[fidx] +=
					factF * (((*swa.edge_f2_)[ieidx] + (*swa.edge_s2R_)[ieidx]) * (*swa.edge_normX_)[ieidx] -
							 (*swa.edge_f3_)[ieidx] * (*swa.edge_normY_)[ieidx]);
				(*swa.face_r_)[fidx] +=
					factF * ((*swa.edge_f3_)[ieidx] * (*swa.edge_normX_)[ieidx] +
							 ((*swa.edge_f2_)[ieidx] + (*swa.edge_s2R_)[ieidx]) * (*swa.edge_normY_)[ieidx]);
			}
			return true;
		});
		return true;
	});

	parallel_foreach_cell(m, [&](Face f) -> bool {
		uint32 fidx = index_of(m, f);

		// friction
		if (swc.friction_ != 0)
		{
			Scalar qx = (*swa.face_q_)[fidx] * cos(swc.alphaK_) + (*swa.face_r_)[fidx] * sin(swc.alphaK_);
			Scalar qy = -(*swa.face_q_)[fidx] * sin(swc.alphaK_) + (*swa.face_r_)[fidx] * cos(swc.alphaK_);
			if ((*swa.face_h_)[fidx] > swc.hmin_)
			{
				qx =
					qx *
					exp(-(9.81 * sqrt(qx * qx + qy * qy) /
						  (std::max(swc.kx_ * swc.kx_, swc.small_ * swc.small_) * pow((*swa.face_h_)[fidx], 7. / 3.))) *
						swc.dt_);
				qy =
					qy *
					exp(-(9.81 * sqrt(qx * qx + qy * qy) /
						  (std::max(swc.ky_ * swc.ky_, swc.small_ * swc.small_) * pow((*swa.face_h_)[fidx], 7. / 3.))) *
						swc.dt_);
			}
			else
			{
				qx = 0.;
				qy = 0.;
			}
			(*swa.face_q_)[fidx] = qx * cos(swc.alphaK_) - qy * sin(swc.alphaK_);
			(*swa.face_r_)[fidx] = qx * sin(swc.alphaK_) + qy * cos(swc.alphaK_);
		}

		// optional correction
		// Negative water depth
		if ((*swa.face_h_)[fidx] < 0.)
		{
			(*swa.face_h_)[fidx] = 0.;
			// (*swa.face_h_)[fidx] = swc.hmin_;
			(*swa.face_q_)[fidx] = 0.;
			(*swa.face_r_)[fidx] = 0.;
		}

		// Abnormal large velocity => Correction of q and r to respect Vmax and Frmax
		if ((*swa.face_h_)[fidx] > swc.hmin_)
		{
			Scalar v = sqrt((*swa.face_q_)[fidx] * (*swa.face_q_)[fidx] + (*swa.face_r_)[fidx] * (*swa.face_r_)[fidx]) /
					   std::max((*swa.face_h_)[fidx], swc.small_);
			Scalar c = sqrt(9.81 * std::max((*swa.face_h_)[fidx], swc.small_));
			Scalar Fr = v / c;
			Scalar Fact = std::max({1.0, v / swc.v_max_, Fr / swc.Fr_max_});
			(*swa.face_q_)[fidx] /= Fact;
			(*swa.face_r_)[fidx] /= Fact;
		}
		else // Quasi-zero
		{
			(*swa.face_q_)[fidx] = 0.;
			(*swa.face_r_)[fidx] = 0.;
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

	swc.t_ += swc.dt_;
	// nb_iter_++;

	auto end = std::chrono::high_resolution_clock::now();

	std::chrono::nanoseconds sleep_duration =
		std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::duration<Scalar>(swc.dt_)) -
		std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

	if (sleep_duration > std::chrono::nanoseconds::zero())
		std::this_thread::sleep_for(sleep_duration);
}

} // namespace shallow_water

} // namespace simulation

} // namespace cgogn

#endif // CGOGN_SIMULATION_SHALLOW_WATER_H_

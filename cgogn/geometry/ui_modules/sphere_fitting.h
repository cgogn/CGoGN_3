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

#ifndef CGOGN_MODULE_SPHERE_FITTING_H_
#define CGOGN_MODULE_SPHERE_FITTING_H_

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/geometry/algos/medial_axis.h>
#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/functions/distance.h>

#include <Eigen/Sparse>
#include <libacc/bvh_tree.h>
#include <libacc/kd_tree.h>

#include <GLFW/glfw3.h>

namespace cgogn
{

namespace ui
{

using geometry::Scalar;
using geometry::Vec3;
using geometry::Vec4;

template <typename SURFACE, typename POINTS>
class SphereFitting : public ViewModule
{
	template <typename T>
	using SAttribute = typename mesh_traits<SURFACE>::template Attribute<T>;
	using SVertex = typename mesh_traits<SURFACE>::Vertex;
	using SFace = typename mesh_traits<SURFACE>::Face;

	template <typename T>
	using PAttribute = typename mesh_traits<POINTS>::template Attribute<T>;
	using PVertex = typename mesh_traits<POINTS>::Vertex;

public:
	SphereFitting(const App& app)
		: ViewModule(app, "SphereFitting (" + std::string{mesh_traits<SURFACE>::name} + "," +
							  std::string{mesh_traits<POINTS>::name} + ")")
	{
	}
	~SphereFitting()
	{
	}

	void fit_sphere(SURFACE& s, const SAttribute<Vec3>* surface_vertex_position,
					const SAttribute<Vec3>* surface_vertex_normal, const acc::BVHTree<uint32, Vec3>* surface_bvh,
					const std::vector<SFace>& surface_bvh_faces, const acc::KDTree<3, uint32>* surface_kdt,
					const std::vector<SVertex>& surface_kdt_vertices, CellsSet<SURFACE, SVertex>* surface_vertices_set,
					POINTS& spheres, PAttribute<Vec3>* spheres_position, PAttribute<Scalar>* spheres_radius,
					PAttribute<Vec4>* spheres_color)
	{
		// 1 - linearized least squares fitting of a sphere to a set of points
		// ----------------------------------------------------------------
		Eigen::MatrixXd A(surface_vertices_set->size(), 4);
		Eigen::VectorXd b(surface_vertices_set->size());
		uint32 idx = 0;
		surface_vertices_set->foreach_cell([&](SVertex v) {
			const Vec3& p = value<Vec3>(s, surface_vertex_position, v);
			A.row(idx) = Eigen::Vector4d(-2.0 * p[0], -2.0 * p[1], -2.0 * p[2], 1.0);
			b(idx) = -(p[0] * p[0]) - (p[1] * p[1]) - (p[2] * p[2]);
			++idx;
		});
		Eigen::LDLT<Eigen::MatrixXd> solver(A.transpose() * A);
		Eigen::MatrixXd s1 = solver.solve(A.transpose() * b);
		s1(3) = std::sqrt(s1(0) * s1(0) + s1(1) * s1(1) + s1(2) * s1(2) - s1(3));

		Vec3 s1c = Vec3(s1(0), s1(1), s1(2));
		Scalar s1r = s1(3);

		Scalar error1 = 0.0;
		surface_vertices_set->foreach_cell([&](SVertex v) {
			Scalar d = (value<Vec3>(s, surface_vertex_position, v) - s1c).norm() - s1r;
			error1 += d * d;
		});
		std::cout << "error1: " << error1 << std::endl;

		// 2 - non-linear optimization of the fitted sphere
		// --------------------------------------------
		Eigen::MatrixXd J(surface_vertices_set->size(), 4);
		Eigen::VectorXd r(surface_vertices_set->size());
		Eigen::VectorXd s2(4);
		s2 << s1(0), s1(1), s1(2), s1(3);
		for (uint32 i = 0; i < 20; ++i)
		{
			idx = 0;
			surface_vertices_set->foreach_cell([&](SVertex v) {
				const Vec3& p = value<Vec3>(s, surface_vertex_position, v);
				Vec3 d = p - Vec3(s2(0), s2(1), s2(2));
				Scalar l = d.norm();
				J.row(idx) = Eigen::Vector4d(-(d[0] / l), -(d[1] / l), -(d[2] / l), -1.0);
				r(idx) = -(l - s2(3));
				++idx;
			});
			Eigen::LDLT<Eigen::MatrixXd> solver(J.transpose() * J);
			s2 += solver.solve(J.transpose() * r);
		}

		Vec3 s2c = Vec3(s2(0), s2(1), s2(2));
		Scalar s2r = s2(3);

		Scalar error2 = 0.0;
		surface_vertices_set->foreach_cell([&](SVertex v) {
			Scalar d = (value<Vec3>(s, surface_vertex_position, v) - s2c).norm() - s2r;
			error2 += d * d;
		});
		std::cout << "error2: " << error2 << std::endl;

		// 3 - medial axis projection using SBC of closest vertex on the surface (with its normal as search axis)
		// --------------------------------------------------------------------------------------------------
		std::pair<uint32, Scalar> k_res;
		surface_kdt->find_nn(s2c, &k_res);
		SVertex closest_vertex = surface_kdt_vertices[k_res.first];
		auto [s3c, s3r, q3] =
			geometry::shrinking_ball_center(s, value<Vec3>(s, surface_vertex_position, closest_vertex),
											value<Vec3>(s, surface_vertex_normal, closest_vertex),
											surface_vertex_position, surface_bvh, surface_bvh_faces, surface_kdt);

		Scalar error3 = 0.0;
		surface_vertices_set->foreach_cell([&](SVertex v) {
			Scalar d = (value<Vec3>(s, surface_vertex_position, v) - s3c).norm() - s3r;
			error3 += d * d;
		});
		std::cout << "error3: " << error3 << std::endl;

		// 4 - same as 3, but use the closest vertex direction as SBC search axis (instead of normal)
		// ------------------------------------------------------------------------------------------
		Vec3 closest_vertex_position = value<Vec3>(s, surface_vertex_position, closest_vertex);
		Vec3 closest_vertex_dir = (closest_vertex_position - s2c).normalized();
		auto [s4c, s4r, q4] =
			geometry::shrinking_ball_center(s, closest_vertex_position, closest_vertex_dir, surface_vertex_position,
											surface_bvh, surface_bvh_faces, surface_kdt);

		Scalar error4 = 0.0;
		surface_vertices_set->foreach_cell([&](SVertex v) {
			Scalar d = (value<Vec3>(s, surface_vertex_position, v) - s4c).norm() - s4r;
			error4 += d * d;
		});
		std::cout << "error4: " << error4 << std::endl;

		// 5 - same as 4, but use the closest point on the surface instead of the closest vertex
		// -------------------------------------------------------------------------------------
		std::pair<uint32, Vec3> bvh_res;
		surface_bvh->closest_point(s2c, &bvh_res);
		Vec3 closest_point_position = bvh_res.second;
		Vec3 closest_point_dir = (closest_point_position - s2c).normalized();
		auto [s5c, s5r, q5] =
			geometry::shrinking_ball_center(s, closest_point_position, closest_point_dir, surface_vertex_position,
											surface_bvh, surface_bvh_faces, surface_kdt);

		Scalar error5 = 0.0;
		surface_vertices_set->foreach_cell([&](SVertex v) {
			Scalar d = (value<Vec3>(s, surface_vertex_position, v) - s5c).norm() - s5r;
			error5 += d * d;
		});
		std::cout << "error5: " << error5 << std::endl;

		// medial axis projection using closest medial axis sample
		// -------------------------------------------------------

		// add spheres
		// -----------
		PVertex v5 = add_vertex(spheres);
		value<Vec3>(spheres, spheres_position, v5) = s5c;
		value<Scalar>(spheres, spheres_radius, v5) = s5r;
		value<Vec4>(spheres, spheres_color, v5) = Vec4(1.0, 1.0, 1.0, 0.5);

		// PVertex v4 = add_vertex(spheres);
		// value<Vec3>(spheres, spheres_position, v4) = s4c;
		// value<Scalar>(spheres, spheres_radius, v4) = s4r;
		// value<Vec4>(spheres, spheres_color, v4) = Vec4(0.0, 0.0, 1.0, 0.4);

		// PVertex v3 = add_vertex(spheres);
		// value<Vec3>(spheres, spheres_position, v3) = s3c;
		// value<Scalar>(spheres, spheres_radius, v3) = s3r;
		// value<Vec4>(spheres, spheres_color, v3) = Vec4(1.0, 0.0, 0.0, 0.4);

		PVertex v2 = add_vertex(spheres);
		value<Vec3>(spheres, spheres_position, v2) = s2c;
		value<Scalar>(spheres, spheres_radius, v2) = s2r;
		value<Vec4>(spheres, spheres_color, v2) = Vec4(0.0, 1.0, 0.0, 0.2);

		// PVertex v1 = add_vertex(spheres);
		// value<Vec3>(spheres, spheres_position, v1) = s1c;
		// value<Scalar>(spheres, spheres_radius, v1) = s1r;
		// value<Vec4>(spheres, spheres_color, v1) = Vec4(1.0, 0.0, 0.0, 0.2);

		// update rendering data
		// ---------------------
		points_provider_->emit_connectivity_changed(spheres);
		points_provider_->emit_attribute_changed(spheres, spheres_position);
		points_provider_->emit_attribute_changed(spheres, spheres_radius);
		points_provider_->emit_attribute_changed(spheres, spheres_color);
	}

	void set_surface(SURFACE& s)
	{
		surface_ = &s;
		surface_vertex_position_.reset();
		surface_vertex_normal_.reset();
		if (surface_bvh_)
		{
			delete surface_bvh_;
			surface_bvh_ = nullptr;
			surface_bvh_faces_.clear();
		}
		if (surface_kdt_)
		{
			delete surface_kdt_;
			surface_kdt_ = nullptr;
			surface_kdt_vertices_.clear();
		}
		// medial_axis_samples_position_.reset();
		// medial_axis_samples_radius_.reset();
		// medial_axis_samples_closest_points_.reset();
		surface_vertices_set_ = nullptr;
	}
	void set_surface_vertex_position(const std::shared_ptr<SAttribute<Vec3>>& surface_vertex_position)
	{
		surface_vertex_position_ = surface_vertex_position;
	}
	void set_surface_vertex_normal(const std::shared_ptr<SAttribute<Vec3>>& surface_vertex_normal)
	{
		surface_vertex_normal_ = surface_vertex_normal;
	}
	void set_surface_vertices_set(CellsSet<SURFACE, SVertex>* surface_vertices_set)
	{
		surface_vertices_set_ = surface_vertices_set;
	}

	void init_surface_data(SURFACE& s, SAttribute<Vec3>* vertex_position, SAttribute<Vec3>* vertex_normal)
	{
		uint32 nb_vertices = nb_cells<SVertex>(s);
		uint32 nb_faces = nb_cells<SFace>(s);

		auto bvh_vertex_index = add_attribute<uint32, SVertex>(s, "__bvh_vertex_index");

		surface_kdt_vertices_.clear();
		surface_kdt_vertices_.reserve(nb_vertices);
		std::vector<Vec3> vertex_position_vector;
		vertex_position_vector.reserve(nb_vertices);
		uint32 idx = 0;
		foreach_cell(s, [&](SVertex v) -> bool {
			surface_kdt_vertices_.push_back(v);
			value<uint32>(s, bvh_vertex_index, v) = idx++;
			vertex_position_vector.push_back(value<Vec3>(s, vertex_position, v));
			return true;
		});

		surface_bvh_faces_.clear();
		surface_bvh_faces_.reserve(nb_faces);
		std::vector<uint32> face_vertex_indices;
		face_vertex_indices.reserve(nb_faces * 3);
		foreach_cell(s, [&](SFace f) -> bool {
			surface_bvh_faces_.push_back(f);
			foreach_incident_vertex(s, f, [&](SVertex v) -> bool {
				face_vertex_indices.push_back(value<uint32>(s, bvh_vertex_index, v));
				return true;
			});
			return true;
		});

		surface_bvh_ = new acc::BVHTree<uint32, Vec3>(face_vertex_indices, vertex_position_vector);
		surface_kdt_ = new acc::KDTree<3, uint32>(vertex_position_vector);

		remove_attribute<SVertex>(s, bvh_vertex_index);
	}

	void set_spheres(POINTS& p)
	{
		spheres_ = &p;
		spheres_position_.reset();
		spheres_radius_.reset();
		spheres_color_.reset();
	}
	void set_spheres_position(const std::shared_ptr<PAttribute<Vec3>>& spheres_position)
	{
		spheres_position_ = spheres_position;
	}
	void set_spheres_radius(const std::shared_ptr<PAttribute<Scalar>>& spheres_radius)
	{
		spheres_radius_ = spheres_radius;
	}
	void set_spheres_color(const std::shared_ptr<PAttribute<Vec4>>& spheres_color)
	{
		spheres_color_ = spheres_color;
	}

protected:
	void init() override
	{
		surface_provider_ = static_cast<ui::MeshProvider<SURFACE>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<SURFACE>::name} + ")"));

		points_provider_ = static_cast<ui::MeshProvider<POINTS>*>(
			app_.module("MeshProvider (" + std::string{mesh_traits<POINTS>::name} + ")"));
	}

	void key_press_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_I)
		{
			int32 x = view->mouse_x();
			int32 y = view->mouse_y();

			rendering::GLVec3d near = view->unproject(x, y, 0.0);
			rendering::GLVec3d far_d = view->unproject(x, y, 1.0);
			Vec3 A{near.x(), near.y(), near.z()};
			Vec3 B{far_d.x(), far_d.y(), far_d.z()};

			Vec3 picked_sphere_center;
			foreach_cell(*spheres_, [&](PVertex v) -> bool {
				if (!picked_sphere_.is_valid())
				{
					picked_sphere_ = v;
					picked_sphere_center = value<Vec3>(*spheres_, spheres_position_, picked_sphere_);
					return true;
				}
				const Vec3& sp = value<Vec3>(*spheres_, spheres_position_, v);
				// Scalar sr = value<Scalar>(*spheres, spheres_radius_, v);
				if (geometry::squared_distance_line_point(A, B, sp) <
					geometry::squared_distance_line_point(A, B, picked_sphere_center))
				{
					picked_sphere_ = v;
					picked_sphere_center = sp;
				}
				return true;
			});

			sphere_info_popup_ = true;
		}
	}

	void key_release_event(View* view, int32 key_code) override
	{
		if (key_code == GLFW_KEY_I)
			sphere_info_popup_ = false;
	}

	void popups() override
	{
		if (sphere_info_popup_)
			ImGui::OpenPopup("Sphere Info");

		if (ImGui::BeginPopup("Sphere Info"))
		{
			if (picked_sphere_.is_valid())
			{
				ImGui::Text("Picked sphere:");
				const Vec3& sp = value<Vec3>(*spheres_, spheres_position_, picked_sphere_);
				ImGui::Text("Center: (%f, %f, %f)", sp[0], sp[1], sp[2]);
				ImGui::Text("Radius: %f", value<Scalar>(*spheres_, spheres_radius_, picked_sphere_));
			}
			else
			{
				ImGui::Text("No sphere picked");
			}
			ImGui::EndPopup();
		}
	}

	void left_panel() override
	{
		imgui_mesh_selector(surface_provider_, surface_, "Surface", [&](SURFACE& s) { set_surface(s); });

		if (surface_)
		{
			MeshData<SURFACE>& md = surface_provider_->mesh_data(*surface_);

			imgui_combo_attribute<SVertex, Vec3>(
				*surface_, surface_vertex_position_, "Position",
				[&](const std::shared_ptr<SAttribute<Vec3>>& attribute) { set_surface_vertex_position(attribute); });

			imgui_combo_attribute<SVertex, Vec3>(
				*surface_, surface_vertex_normal_, "Normal",
				[&](const std::shared_ptr<SAttribute<Vec3>>& attribute) { set_surface_vertex_normal(attribute); });

			imgui_combo_cells_set(md, surface_vertices_set_, "Vertices",
								  [&](CellsSet<SURFACE, SVertex>* cs) { set_surface_vertices_set(cs); });
		}

		imgui_mesh_selector(points_provider_, spheres_, "Spheres", [&](POINTS& p) { set_spheres(p); });

		if (spheres_)
		{
			MeshData<POINTS>& md = points_provider_->mesh_data(*spheres_);

			imgui_combo_attribute<PVertex, Vec3>(
				*spheres_, spheres_position_, "Position",
				[&](const std::shared_ptr<PAttribute<Vec3>>& attribute) { set_spheres_position(attribute); });

			imgui_combo_attribute<PVertex, Scalar>(
				*spheres_, spheres_radius_, "Radius",
				[&](const std::shared_ptr<PAttribute<Scalar>>& attribute) { set_spheres_radius(attribute); });

			imgui_combo_attribute<PVertex, Vec4>(
				*spheres_, spheres_color_, "Color",
				[&](const std::shared_ptr<PAttribute<Vec4>>& attribute) { set_spheres_color(attribute); });
		}

		ImGui::Separator();

		if (surface_ && surface_vertex_position_ && surface_vertex_normal_)
		{
			if (ImGui::Button("Init surface data"))
				init_surface_data(*surface_, surface_vertex_position_.get(), surface_vertex_normal_.get());
		}

		if (surface_ && surface_vertex_position_ && surface_vertex_normal_ && surface_bvh_ && surface_kdt_ &&
			surface_vertices_set_ && spheres_ && spheres_position_ && spheres_radius_ && spheres_color_)
		{
			if (ImGui::Button("Fit sphere"))
			{
				fit_sphere(*surface_, surface_vertex_position_.get(), surface_vertex_normal_.get(), surface_bvh_,
						   surface_bvh_faces_, surface_kdt_, surface_kdt_vertices_, surface_vertices_set_, *spheres_,
						   spheres_position_.get(), spheres_radius_.get(), spheres_color_.get());
			}
		}

		ImGui::Separator();

		if (picked_sphere_.is_valid())
		{
			ImGui::Text("Picked sphere:");
			const Vec3& sp = value<Vec3>(*spheres_, spheres_position_, picked_sphere_);
			ImGui::Text("Center: (%f, %f, %f)", sp[0], sp[1], sp[2]);
			ImGui::Text("Radius: %f", value<Scalar>(*spheres_, spheres_radius_, picked_sphere_));
		}
		else
		{
			ImGui::Text("No sphere picked");
		}
	}

private:
	SURFACE* surface_ = nullptr;
	std::shared_ptr<SAttribute<Vec3>> surface_vertex_position_ = nullptr;
	std::shared_ptr<SAttribute<Vec3>> surface_vertex_normal_ = nullptr;
	acc::BVHTree<uint32, Vec3>* surface_bvh_ = nullptr;
	std::vector<SFace> surface_bvh_faces_;
	acc::KDTree<3, uint32>* surface_kdt_ = nullptr;
	std::vector<SVertex> surface_kdt_vertices_;
	// std::shared_ptr<SAttribute<Vec3>> medial_axis_samples_position_ = nullptr;
	// std::shared_ptr<SAttribute<Scalar>> medial_axis_samples_radius_ = nullptr;
	// std::shared_ptr<SAttribute<std::pair<Vec3, Vec3>>> medial_axis_samples_closest_points_ = nullptr;
	CellsSet<SURFACE, SVertex>* surface_vertices_set_ = nullptr;

	POINTS* spheres_ = nullptr;
	std::shared_ptr<PAttribute<Vec3>> spheres_position_ = nullptr;
	std::shared_ptr<PAttribute<Scalar>> spheres_radius_ = nullptr;
	std::shared_ptr<PAttribute<Vec4>> spheres_color_ = nullptr;

	MeshProvider<SURFACE>* surface_provider_ = nullptr;
	MeshProvider<POINTS>* points_provider_ = nullptr;

	bool sphere_info_popup_ = false;
	PVertex picked_sphere_;
};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_SPHERE_FITTING_H_

/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
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

#ifndef CGOGN_RENDERING_SKEL_SHAPE_DRAWER_H_
#define CGOGN_RENDERING_SKEL_SHAPE_DRAWER_H_

#include <cgogn/rendering/cgogn_rendering_export.h>
#include <cgogn/rendering/shader_program.h>
#include <cgogn/rendering/vbo_update.h>

#include <cgogn/geometry/types/vector_traits.h>
#include <Eigen/Dense>

namespace cgogn
{

namespace rendering
{
/*
ajouter:
void ShaderParam::set_vbos_inst(const std::vector<VBO*>& vbos)
void VBO::associate(GLuint attrib, int32 stride = 0, uint32 first = 0, uint32 divisor = 0 )
*/	
class ShaderParamSkelShape : public ShaderParam
{
public:

	//subdiv
	int nbs_;
	//color RBGA 
	GLColor color_;

	GLVec3 light_position_;

	std::unique_ptr<rendering::VAO> inst_vao_;

	ShaderParamSkelShape(ShaderProgram* prg);

	virtual void draw_inst(const GLMat4& projection, const GLMat4& view, uint32 nbi) = 0;

	inline void set_subdiv(int32 nb)
	{
		nbs_ = nb;
	}

	inline void set_light_pos(const GLVec3& lp)
	{
		light_position_ = lp;
	}

	inline void set_color(const GLColor& col)
	{
		color_ = col;
	}

};


DECLARE_SHADER_CLASS(SkelSphere, false, CGOGN_STR(SkelSphere))

class CGOGN_RENDERING_EXPORT ShaderParamSkelSphere : public ShaderParamSkelShape
{
	void set_uniforms() override;

public:

	using ShaderType = ShaderSkelSphere;

	inline ShaderParamSkelSphere(ShaderType* sh):
	ShaderParamSkelShape(sh) 
	{}

	inline ~ShaderParamSkelSphere() override {}

	void draw_inst(const GLMat4& projection, const GLMat4& view, uint32 nbi) override;
};

class SkelSphereDrawer
{
	VBO vbo_;
public:
	std::unique_ptr<ShaderSkelSphere::Param> param_;
	SkelSphereDrawer();
	void set_vertices_spheres(const std::vector<GLVec4>& P1);
	void set_subdiv(uint32 sub);
	void draw(const GLMat4& projection, const GLMat4& view);
};



DECLARE_SHADER_CLASS(SkelCone, false, CGOGN_STR(SkelCone))

class CGOGN_RENDERING_EXPORT ShaderParamSkelCone : public ShaderParamSkelShape
{

	void set_uniforms() override;

public:

	using ShaderType = ShaderSkelCone;

	inline ShaderParamSkelCone(ShaderType* sh):
	ShaderParamSkelShape(sh) 
	{}

	inline ~ShaderParamSkelCone() override {}

	void draw_inst(const GLMat4& projection, const GLMat4& view, uint32 nbi) override;
};

class SkelConeDrawer
{
	VBO vbo1_;
	VBO vbo2_;
	VBO vbo3_;
	VBO vbo4_;

	void set_real_cones(const std::vector<GLVec4>& P1, const std::vector<GLVec4>& P2, const std::vector<GLVec3>& N1,
						const std::vector<GLVec3>& N2);
	void compute_skel_cone(const GLVec4& A, const GLVec4& B, GLVec4& P1, GLVec4& P2, GLVec3& N1, GLVec3& N2);

public:
	std::unique_ptr<ShaderSkelCone::Param> param_;

	SkelConeDrawer();

	void set_edges_cones(const std::vector<GLVec4>& edges);

	void draw(const GLMat4& projection, const GLMat4& view);

};


DECLARE_SHADER_CLASS(SkelTri, false, CGOGN_STR(SkelTri))

class CGOGN_RENDERING_EXPORT ShaderParamSkelTri : public ShaderParamSkelShape
{
	void set_uniforms() override;
public:

	using ShaderType = ShaderSkelTri;

	inline ShaderParamSkelTri(ShaderType* sh):
	ShaderParamSkelShape(sh) 
	{}

	inline ~ShaderParamSkelTri() override {}

	void draw_inst(const GLMat4& projection, const GLMat4& view, uint32 nbi) override;
};

class SkelTriDrawer
{
	VBO vbo_p_;
	VBO vbo_n_;

	void compute_CDisk(const GLVec4& A, const GLVec4& B, GLVec3& P1, GLVec3& P2);
	void interPPS(const GLVec4& Pa, const GLVec4& Pb,
					const GLVec4& Sph, GLVec3& I1, GLVec3& I2);
	void compute_tri2(const GLVec4& P1, const GLVec4& P2, const GLVec4& P3, GLVec3* I);

public:
	std::unique_ptr<ShaderSkelTri::Param> param_;
	SkelTriDrawer();
	void set_tris_faces(const std::vector<GLVec4>& pts);
	void draw(const GLMat4& projection, const GLMat4& view);
};


/*
  must be created with correct OGL context !
  usage:
  1 fill the buffers with add_vertex add_edge and add tri
    (parameter can be any type of VEC which offer [] operator and good number of components.
  2 call update
  3 draw
*/
class SkelShapeDrawer
{
	std::vector<GLVec4> vertices_;
	std::vector<GLVec4> edges_;
	std::vector<GLVec4> triangles_;

	SkelSphereDrawer sphere_drawer_;
	SkelConeDrawer cone_drawer_;
	SkelTriDrawer tri_drawer_;

public:
	inline SkelShapeDrawer()
	{}


	template<typename VEC4>
	void add_vertex(const VEC4& v)
	{
		vertices_.emplace_back(float(v[0]), float(v[1]), float(v[2]), float(v[3]));
	}

	template <typename VEC3, typename SCAL>
	inline void add_vertex(const VEC3& v, SCAL r)
	{
		vertices_.emplace_back(float(v[0]), float(v[1]), float(v[2]), float(r));
	}

	template <typename VEC4_1, typename VEC4_2>
	void add_edge(const VEC4_1& e1, const VEC4_2& e2)
	{
		edges_.emplace_back(float(e1[0]), float(e1[1]), float(e1[2]), float(e1[3]));
		edges_.emplace_back(float(e2[0]), float(e2[1]), float(e2[2]), float(e2[3]));
	}

	template <typename VEC3_1, typename VEC3_2, typename SCAL_1, typename SCAL_2>
	void add_edge(const VEC3_1& e1, SCAL_1 r1, const VEC3_2& e2, SCAL_2 r2)
	{
		edges_.emplace_back(float(e1[0]), float(e1[1]), float(e1[2]), float(r1));
		edges_.emplace_back(float(e2[0]), float(e2[1]), float(e2[2]), float(r2));
	}

	template <typename VEC4_1, typename VEC4_2, typename VEC4_3>
	void add_triangle(const VEC4_1& t1, const VEC4_2& t2, const VEC4_3& t3)
	{
		triangles_.emplace_back(float(t1[0]), float(t1[1]), float(t1[2]), float(t1[3]));
		triangles_.emplace_back(float(t2[0]), float(t2[1]), float(t2[2]), float(t2[3]));
		triangles_.emplace_back(float(t3[0]), float(t3[1]), float(t3[2]), float(t3[3]));
	}

	template <typename VEC3_1, typename VEC3_2, typename VEC3_3, typename SCAL_1, typename SCAL_2, typename SCAL_3>
	void add_triangle(const VEC3_1& t1, SCAL_1 r1, const VEC3_2& t2, SCAL_2 r2, const VEC3_3& t3, SCAL_3 r3)
	{
		triangles_.emplace_back(float(t1[0]), float(t1[1]), float(t1[2]), float(r1));
		triangles_.emplace_back(float(t2[0]), float(t2[1]), float(t2[2]), float(r2));
		triangles_.emplace_back(float(t3[0]), float(t3[1]), float(t3[2]), float(r3));
	}

	/*
	 * eed to be called after filling the buffers (and before calling draw)
	 */
	inline void update()
	{
		sphere_drawer_.set_vertices_spheres(vertices_);
		cone_drawer_.set_edges_cones(edges_);
		tri_drawer_.set_tris_faces(triangles_);
	}

	inline void draw(const GLMat4& proj_matrix, const GLMat4& view_matrix)
	{
		sphere_drawer_.draw(proj_matrix, view_matrix);
		cone_drawer_.draw(proj_matrix, view_matrix);
		tri_drawer_.draw(proj_matrix, view_matrix);
	}

	inline void set_subdiv(uint32 nbs)
	{
		sphere_drawer_.set_subdiv(nbs);
		cone_drawer_.param_->set_subdiv(nbs);
	}

	inline void set_color(const GLColor col)
	{
		sphere_drawer_.param_->set_color(col);
		cone_drawer_.param_->set_color(col);
		tri_drawer_.param_->set_color(col);
	}
};


template <typename VEC4,typename VEC3, typename SCALAR>
class SkeletonSampler
{
public:

	void add_vertex(const VEC4& v)
	{
		//	vertices_.emplace_back((VEC4(v[0], v[1], v[2], v[3]);
		vertices_push_back(v);
	}

	void add_vertex(const VEC3& v, SCALAR r)
	{
	//		vertices_.emplace_back(stdd::make_pair(VEC(double(v[0]), double(v[1]), double(v[2])), SCALAR(r)));
		vertices_.emplace_back(SCALAR(v[0]),SCALAR(v[1]),SCALAR(v[2]),SCALAR(r));
	}


	void add_edge(const VEC3& A, SCALAR Ra, const VEC3& B, SCALAR Rb)
	{
		// copmpute real cone center and radius
		VEC3 AB = B - A;

		SCALAR d = AB.norm();
		SCALAR k = (Rb - Ra) / (d * d);
		SCALAR kk = std::sqrt(1.0 - (Rb - Ra) * k);

		VEC3 Ea = A - AB * (k * Ra);
		VEC3 Eb = B - AB * (k * Rb);
		edges_.emplace_back(Ea[0], Ea[1], Ea[2], Ra * kk);
		edges_.emplace_back(Eb[0], Eb[1], Eb[2], Rb * kk);
	}

	inline void add_edge(const VEC4& A, const VEC4& B)
	{
		add_edge(A.template topRows<3>(), A[3], B.template topRows<3>(), B[3]);
	}

	void add_triangle(const VEC3& A, SCALAR Ra, const VEC3& B, SCALAR Rb, const VEC3& C, SCALAR Rc)
	{
	}

	inline void add_triangle(const VEC4& A, const VEC4& B, const VEC4& C)
	{
		auto nb = triangles_.size();
		triangles_.resize(nb+5);
		compute_prism(A,B,C,&(triangles_[nb]));
	}

	inline SCALAR eval_skeketon(const VEC3& P)
	{
		auto it = vertices_.begin();
		auto dist = evalSphereSDF(P, it++);
		while (it != vertices_.end())
			dist = std::min(dist,evalSphereSDF(P,++it));
		
		it = edges_.begin();
		while (it != edges_.end())
			dist = std::min(dist,evalConeSDF(P,it++++));

		it = triangles_.begin();
		while (it != triangles_.end())
		{
			dist = std::min(dist,evalConeSDF(P,it));
			it += 5;
		}
	}
	
protected:
	std::vector<VEC4> vertices_;	 // one by one
	std::vector<VEC4> edges_;	 // by pair
	std::vector<VEC4> triangles_; // by 5

	static inline SCALAR evalSphereSDF( const VEC3& P, typename std::vector<VEC4>::const_iterator it)
	{
		return (it->template topRows<3>() - P).norm() - it->w();
	}

	static inline SCALAR evalConeSDF(
		const VEC3& P, typename std::vector<VEC4>::const_iterator it)
	{
		const VEC3& P_a = it->template topRows<3>();
		SCALAR ra = it->w();
		++it;
		const VEC3& Pb = it->first;
		SCALAR rb = it->second;

		SCALAR rba = rb - ra;
		VEC3 ba = Pb - P_a;
		SCALAR baba = dot(ba, ba);
		VEC3 pa = P - P_a;
		SCALAR papa = dot(pa, pa);
		SCALAR paba = dot(pa, ba) / baba;
		SCALAR x = std::sqrt(papa - paba * paba * baba);
		SCALAR cax = std::max(SCALAR(0), x - ((paba < SCALAR(0.5)) ? ra : rb));
		SCALAR cay = std::abs(paba - 0.5) - 0.5;
		SCALAR k = rba * rba + baba;
		SCALAR f = (rba * (x - ra) + paba * baba) / k;
		if (f < SCALAR(0))
			f = SCALAR(0);
		if (f > SCALAR(1))
			f = SCALAR(1);
		SCALAR cbx = x - ra - f * rba;
		SCALAR cby = paba - f;
		SCALAR s = (cbx < SCALAR(0) && cay < SCALAR(0)) ? SCALAR(-1) : SCALAR(1);
		return s * std::sqrt(std::min(cax * cax + cay * cay * baba, cbx * cbx + cby * cby * baba));
	}

	static inline SCALAR evalPlaneSDF(const VEC3& P, typename std::vector<VEC4>::const_iterator it)
	{
		return P.dot(it->template topRows<3>()) + it->w();
	}

	static inline SCALAR evalPrismTriSDF(const VEC3& P, typename std::vector<VEC4>::const_iterator it)
	{
		SCALAR dist = evalPlaneSDF(P,it);
		for(int i=1;i<6;++i)
			dist = std::max(dist,evalPlaneSDF(P,++it));
		return dist;
	}

	static inline void compute_CenterDisk(const VEC4& sphA, const VEC4& sphB, VEC3& centerA, VEC3& centerB)
	{
		auto AB = sphB - sphA;
		const VEC3& AB3 = AB.template topRows<3>();
		auto d = AB3.norm();
		auto k = AB[3] / (d * d);
		centerA = sphA.template topRows<3>() - AB3 * (k * sphA[3]);
		centerB = sphB.template topRows<3>() - AB3 * (k * sphB[3]);
	}

	static inline void interPPS(const VEC4& planeA, const VEC4& planeB, const VEC4& Sph, VEC3& I1, VEC3& I2)
	{
		const VEC3& Na = planeA.template topRows<3>();
		const VEC3& Nb = planeB.template topRows<3>();

		VEC3 U = Na.cross(Nb).normalized();

		SCALAR dp = Na.dot(Nb);
		SCALAR c1 = (-planeA[3] + planeB[3] * dp);
		SCALAR c2 = (-planeB[3] + planeA[3] * dp);
		VEC3 O = (c1 * Na + c2 * Nb) / (1.0f - dp * dp);

		VEC3 CO = O - Sph.template topRows<3>();
		dp = U.dot(CO);
		SCALAR delta = dp * dp - (CO.dot(CO) - Sph[3] * Sph[3]);
		SCALAR k1 = -dp + std::sqrt(delta);
		SCALAR k2 = -dp - std::sqrt(delta);

		I1 = O + k1 * U;
		I2 = O + k2 * U;
	};

	void compute_prism(const VEC4& Sph1, const VEC4& Sph2, const VEC4& Sph3, VEC4* Planes)
	{
		std::array<VEC3, 6> centers;

		compute_CDisk(Sph1, Sph2, centers[1], centers[2]);
		VEC3 N1 = (Sph2.template topRows<3>() - Sph1.template topRows<3>()).normalized();
		compute_CDisk(Sph2, Sph3, centers[3], centers[4]);
		VEC3 N2 = (Sph3.template topRows<3>() - Sph2.template topRows<3>()).normalized();
		compute_CDisk(Sph3, Sph1, centers[5], centers[0]);
		VEC3 N3 = (Sph1.template topRows<3>() - Sph3.template topRows<3>()).normalized();

		std::array<VEC3, 6> I;

		interPPS({N1[0], N1[1], N1[2], -N1.dot(centers[1])}, {-N3[0], -N3[1], -N3[2], N3.dot(centers[0])}, Sph1, I[0],
				 I[3]);
		interPPS({N2[0], N2[1], N2[2], -N2.dot(centers[3])}, {-N1[0], -N1[1], -N1[2], N1.dot(centers[2])}, Sph2, I[1],
				 I[5]);
		interPPS({N3[0], N3[1], N3[2], -N3.dot(centers[5])}, {-N2[0], -N2[1], -N2[2], N2.dot(centers[4])}, Sph3, I[2],
				 I[4]);


		Planes[0].template topRows<3>() = (I[1] - I[0]).cross(I[2] - I[0]).normalized();
		Planes[0].w() = Planes[0].template topRows<3>().dot(I[0]);

		Planes[1].template topRows<3>() = (I[4] - I[3]).cross(I[5] - I[3]).normalized();
		Planes[1].w() = Planes[1].template topRows<3>().dot(I[3]);

		Planes[2].template topRows<3>() = (I[5] - I[3]).cross(I[0] - I[3]).normalized();
		Planes[2].w() = Planes[2].template topRows<3>().dot(I[3]);

		Planes[3].template topRows<3>() = (I[4] - I[5]).cross(I[1] - I[5]).normalized();
		Planes[3].w() = Planes[3].template topRows<3>().dot(I[5]);

		Planes[4].template topRows<3>() = (I[3] - I[4]).cross(I[2] - I[4]).normalized();
		Planes[4].w() = Planes[4].template topRows<3>().dot(I[4]);

	}
};




template <typename VEC>
class SkelSDF
{
//public:
//	using SCALAR = typename geometry::vector_traits<VEC>::Scalar
//	inline SkelSDF()
//	{}	using SCALAR = typename VEC::SCALAR;
//
//protected
//	std::vector<std::pair<VEC4> vertices_; // one by one
//	std::vector<std::pair<VEC4> edges_; // by pair
//	std::vector<std::pair<VEC4> triangles_prism_; // by 5
//
//public:
//
//	template<typename VEC4>
//	void add_vertex(const VEC4& v)
//	{
////		vertices_.emplace_back(stdd::make_pair(VEC(double(v[0]), double(v[1]), double(v[2])), SCALAR(v[3])));
//		vertices_.emplace_back((v[0], v[1], v[2]), geometry::vector_traits<VEC>::Scalar(v[3])));
//	}
//
//	template <typename VEC3>
//	inline void add_vertex(const VEC3& v, geometry::vector_traits<VEC3>::Scalar r)
//	{
////		vertices_.emplace_back(stdd::make_pair(VEC(double(v[0]), double(v[1]), double(v[2])), SCALAR(r)));
//		vertices_.emplace_back(stdd::make_pair(VEC (v[0], v[1], v[2]), SCALAR(r)));
//	}


	//template <typename VEC3>
	//void add_edge(const VEC3& A, geometry::vector_traits<VEC3>::Scalar Ra, const VEC3& B,
	//		  geometry::vector_traits<VEC3>::Scalar Rb)
	//{
	//	// copmpute real cone center and radius
	//	auto AB = B - A;

	//	auto d = AB3.norm();
	//	auto k = (Rb-Ra) / (d * d);
	//	auto kk = std::sqrt(1.0 - ((Rb-Ra) * AB[3] / (d * d)));
	//	edges_.emplace_back(stdd::make_pair(A - AB * (k * Ra), Ra * kk);
	//	edges_.emplace_back(stdd::make_pair(B - AB * (k * Rb), Rb * kk);
	//}

	//template <typename VEC4>
	//void add_edge(const VEC4& A, const VEC4& B)
	//{
	//	add_edge(A.template topRows<3>(), A[3], B.template topRows<3>(), B[3]);
	//}

	/*template <typename VEC4, typename VEC3>
	void compute_CenterDisk(const VEC4& sphA, const VEC4& sphB, VEC3& centerA, VEC3& centerB)
	{
		auto AB = sphB - sphA;
		const auto3& AB3 = AB.template topRows<3>();
		auto d = AB3.norm();
		auto k = AB[3] / (d * d);
		centerA = sphA.template topRows<3>() - AB3 * (k * sphA[3]);
		centerB = sphB.template topRows<3>() - AB3 * (k * sphB[3]);
	};*/



	//	auto c1 = (-planeA[3] + planeB[3] * dp);
		//	auto c2 = (-planeB[3] + planeA[3] * dp);
		//	VEC3 O = (c1 * Na + c2 * Nb) / (1.0f - dp * dp);

		//	VEC3 CO = O - Sph.template topRows<3>();
		//	dp = U.dot(CO);
		//	auto delta = dp * dp - (CO.dot(CO) - Sph[3] * Sph[3]);
		//	auto k1 = -dp + std::sqrt(delta);
		//	auto k2 = -dp - std::sqrt(delta);

		//	I1 = O + k1 * U;
		//	I2 = O + k2 * U;
		//};








	template <typename VEC4_1, typename VEC4_2, typename VEC4_3>
	void add_triangle(const VEC4_1& t1, const VEC4_2& t2, const VEC4_3& t3)
	{
// calculer
		// triangles_.emplace_back(stdd::make_pair(VEC(double(t1[0]), double(t1[1]), double(t1[2])), SCALAR(t1[3])));
		// triangles_.emplace_back(stdd::make_pair(VEC(double(t2[0]), double(t2[1]), double(t2[2])), SCALAR(t2[3])));
		// triangles_.emplace_back(stdd::make_pair(VEC(double(t3[0]), double(t3[1]), double(t3[2])), SCALAR(t3[3])));
	}

	template <typename VEC3_1, typename VEC3_2, typename VEC3_3, typename SCAL_1, typename SCAL_2, typename SCAL_3>
	void add_triangle(const VEC3_1& t1, SCAL_1 r1, const VEC3_2& t2, SCAL_2 r2, const VEC3_3& t3, SCAL_3 r3)
	{
// calculer
		// triangles_.emplace_back(stdd::make_pair(VEC(double(t1[0]), double(t1[1]), double(t1[2])), SCALAR(r1)));
		// triangles_.emplace_back(stdd::make_pair(VEC(double(t2[0]), double(t2[1]), double(t2[2])), SCALAR(r2)));
		// triangles_.emplace_back(stdd::make_pair(VEC(double(t3[0]), double(t3[1]), double(t3[2])), SCALAR(r3)));
	}


};



} // namespace rendering

} // namespace cgogn

#endif // CGOGN_RENDERING_SHAPE3_DRAWER_H_





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

#ifndef CGOGN_RENDERING_SKELETON_SAMPLER_H_
#define CGOGN_RENDERING_SKELETON_SAMPLER_H_

#include <cgogn/geometry/types/vector_traits.h>
// #include <Eigen/Dense>

namespace cgogn
{

namespace modeling
{

template <typename VEC4, typename VEC3, typename SCALAR>
class SkeletonSampler
{
public:
	SkeletonSampler()
		: bb_min_(std::numeric_limits<SCALAR>::max(), std::numeric_limits<SCALAR>::max(),
				  std::numeric_limits<SCALAR>::max()),
		  bb_max_(std::numeric_limits<SCALAR>::min(), std::numeric_limits<SCALAR>::min(),
				  std::numeric_limits<SCALAR>::min())
	{
	}

	void add_vertex(const VEC4& v)
	{
		vertices_.push_back(v);

		update_BB(v);
	}

	void add_vertex(const VEC3& v, SCALAR r)
	{
		//		vertices_.emplace_back(stdd::make_pair(VEC(double(v[0]), double(v[1]), double(v[2])), SCALAR(r)));
		vertices_.emplace_back(v[0], v[1], v[2], r);
	}

	void add_edge(const VEC3& A, SCALAR Ra, const VEC3& B, SCALAR Rb)
	{
		// compute real cone center and radius
		VEC3 AB = B - A;

		SCALAR d = AB.norm();
		SCALAR k = (Rb - Ra) / (d * d);
		SCALAR kk = std::sqrt(1.0 - (Rb - Ra) * k);

		VEC3 Ea = A - AB * (k * Ra);
		VEC3 Eb = B - AB * (k * Rb);
		edges_.emplace_back(Ea[0], Ea[1], Ea[2], Ra * kk);
		edges_.emplace_back(Eb[0], Eb[1], Eb[2], Rb * kk);

		update_BB(A, Ra);
		update_BB(B, Rb);
	}

	inline void add_edge(const VEC4& A, const VEC4& B)
	{
		add_edge(A.template topRows<3>(), A[3], B.template topRows<3>(), B[3]);
	}

	inline void add_triangle(const VEC4& A, const VEC4& B, const VEC4& C)
	{
		auto nb = triangles_.size();
		triangles_.resize(nb + 5);
		compute_prism(A, B, C, &(triangles_[nb]));

		update_BB(A);
		update_BB(B);
		update_BB(C);
	}

	inline SCALAR eval_skeleton(const VEC3& P)
	{
		SCALAR dist = std::numeric_limits<SCALAR>::max();

		for (auto it = vertices_.begin(); it != vertices_.end(); ++it)
			dist = std::min(dist, evalSphereSDF(P, it));

		for (auto it = edges_.begin(); it != edges_.end(); it += 2)
			dist = std::min(dist, evalConeSDF(P, it));

		for (auto it = triangles_.begin(); it != triangles_.end(); it += 5)
			dist = std::min(dist, evalPrismTriSDF(P, it));

		return dist;
	}

	// inline void inter_skeleton(const VEC3& P, const VEC3& Du, std::vector<VEC3>& I, SCALAR ds_max)
	// {
	// 	std::cout << "INTER! P=" << P.transpose() << std::endl;
	// 	SCALAR d = eval_skeleton(P);
	// 	SCALAR ds = d;
	// 	VEC3 Q;
	// 	while (ds < ds_max)
	// 	{
	// 		std::cout << "ds:" << ds << " < s_max:" << ds_max << std::endl;
	// 		while ((d > Epsilon) && (ds < ds_max))
	// 		{
	// 			std::cout << "outside " << std::endl;
	// 			std::cout << "d: " << d << std::endl;
	// 			Q = P + ds * Du;
	// 			std::cout << "-> ds: " << ds << " -> Q: " << Q.transpose() << std::endl;
	// 			d = eval_skeleton(Q);
	// 			ds += d;
	// 		}

	// 		SCALAR ds_prec;
	// 		SCALAR d_prec;
	// 		while ((ds < ds_max) && (d <= Epsilon))
	// 		{
	// 			d_prec = d;
	// 			ds_prec = ds;
	// 			Q = P + ds * Du;
	// 			d = eval_skeleton(Q);
	// 			std::cout << "go nearest d=" << d<< std::endl;
	// 			ds += Epsilon;
	// 			if (d <= SCALAR(0))
	// 			{
	// 				std::cout << "ENTER" << d<< std::endl;
	// 				//	I.push_back(P + Du*(ds_prec+ds)/SCALAR(2));
	// 				SCALAR a = (ds * d_prec - ds_prec * d)/(d_prec-d);
	// 				I.push_back(P + Du * a);
	// 				std::cout << "d: " << d_prec << " => " << d << "ds: " << ds_prec << " => " << ds << "=> a = " << a
	// 						  << std::endl;
	// 				std::cout << "PUSH FRONT " << I.back().transpose() << std::endl;
	// 				ds += Epsilon;
	// 				break;
	// 			}
	// 		}

	// 		while ((ds < ds_max) && (-d <= Epsilon))
	// 		{
	// 			d_prec = d;
	// 			ds_prec = ds;
	// 			Q = P + ds * Du;
	// 			d = eval_skeleton(Q);
	// 			ds += Epsilon;
	// 		}

	// 		while ((-d < Epsilon) && (ds < ds_max))
	// 		{
	// 			std::cout << "continue inside " << std::endl;
	// 			std::cout << "d: " << d << std::endl;
	// 			Q = P + ds * Du;
	// 			std::cout << "-> ds: " << ds << " -> Q: " << Q.transpose() << std::endl;
	// 			d = eval_skeleton(Q);
	// 			ds -= d;
	// 		}

	// 		while ((ds < ds_max) && (-d > Epsilon))
	// 		{
	// 			d_prec = d;
	// 			ds_prec = ds;
	// 			Q = P + ds * Du;
	// 			d = eval_skeleton(Q);
	// 			std::cout << "go nearest out d=" << d << std::endl;
	// 			ds += Epsilon;
	// 			if (d >= SCALAR(0))
	// 			{
	// 				std::cout << "LEAVE" << std::endl;

	// 				//					I.push_back(P + Du*(ds_prec+ds)/SCALAR(2));
	// 				SCALAR a = (d * ds_prec - d_prec * ds)/(d-d_prec);
	// 				I.push_back(P + Du * a);
	// 				std::cout << "d: " << d_prec << " => " << d << "ds: " << ds_prec << " => " << ds << "=> a = " << a
	// 						  << std::endl;
	// 				std::cout << "PUSH BACK" << I.back().transpose() << std::endl;
	// 				ds += Epsilon;
	// 				break;
	// 			}
	// 		}

	// 		// d = eval_skeleton(Q);
	// 	}
	// }

	inline void inter_skeleton(const VEC3& P, const VEC3& Du, std::vector<VEC3>& I, SCALAR ds_max)
	{
		// std::cout << "INTER! P=" << P.transpose() << std::endl;
		SCALAR d = eval_skeleton(P);
		SCALAR ds = SCALAR(0);
		SCALAR d_prec = SCALAR(0);;
		SCALAR ds_prec = SCALAR(0);;

		VEC3 Q;
		while (ds < ds_max)
		{
			// std::cout << "ds:" << ds << " < s_max:" << ds_max << std::endl;
			while ((d >= SCALAR(0)) && (ds < ds_max))
			{
				ds_prec = ds;
				ds += std::max(Epsilon,d);
				Q = P + ds * Du;
				d_prec = d;
				d = eval_skeleton(Q);
				// std::cout <<"d="<< d<< " -> ds: " << ds << " -> Q: " << Q.transpose() << std::endl;				
			}
			if (ds >= ds_max)
				break;
			// std::cout << "ENTER OBJ" << d<< std::endl;
			SCALAR a = (ds * d_prec - ds_prec * d)/(d_prec-d);
			I.push_back(P + Du * a);
			// std::cout << "d: " << d_prec << " => " << d << "ds: " << ds_prec << " => " << ds << "=> a = " << a  << std::endl;
			// std::cout << "PUSH FRONT " << I.back().transpose() << std::endl;
//			ds += Epsilon;			

			while ((d < SCALAR(0)) && (ds < ds_max))
			{
				// std::cout << "IN OBJ" << d<< std::endl;
				ds_prec = ds;
				ds += std::max(Epsilon,-d);
				Q = P + ds * Du;
				d_prec = d;
				d = eval_skeleton(Q);
				// std::cout << "-> ds: " << ds << " -> Q: " << Q.transpose() << "  d: "<<d << std::endl;				
			}

			// std::cout << "OUT OBJ" << d<< std::endl;
			a = (ds * d_prec - ds_prec * d)/(d_prec-d);
			I.push_back(P + Du * a);
			// std::cout << "d: " << d_prec << " => " << d << "ds: " << ds_prec << " => " << ds << "=> a = " << a  << std::endl;
			// std::cout << "PUSH BACK " << I.back().transpose() << std::endl;
		}
	}


	// inline void sample(SCALAR step)
	// {
	// 	std::cout << "BB: [" << bb_min_.transpose() << " - " << bb_max_.transpose() << "]" << std::endl;
	// 	Epsilon = step / SCALAR(10);
	// 	samples_.clear();
	// 	samples_.reserve(8192);
	// 	// Z
	// 	for (SCALAR y = bb_min_[1] + step / 2; y < bb_max_[1] - step / 2; y += step)
	// 		for (SCALAR x = bb_min_[0] + step / 2; x < bb_max_[0] - step / 2; x += step)
	// 			inter_skeleton(VEC3{x, y, bb_min_[2]}, VEC3{0, 0, 1}, samples_, bb_max_[2] - bb_min_[2]);

	// 	// Y
	// 	for (SCALAR x = bb_min_[0] + step / 2; x < bb_max_[0] - step / 2; x += step)
	// 		for (SCALAR z = bb_min_[2] + step / 2; z < bb_max_[2] - step / 2; z += step)
	// 			inter_skeleton(VEC3{x, bb_min_[1], z}, VEC3{0, 1, 0}, samples_, bb_max_[1] - bb_min_[1]);

	// 	// X
	// 	for (SCALAR y = bb_min_[1] + step / 2; y < bb_max_[1] - step / 2; y += step)
	// 		for (SCALAR z = bb_min_[2] + step / 2; z < bb_max_[2] - step / 2; z += step)
	// 			inter_skeleton(VEC3{bb_min_[0], y, z}, VEC3{1, 0, 0}, samples_, bb_max_[0] - bb_min_[0]);

	// 	std::cout << samples_.size() << "points genrated" << std::endl;
	// }

	void sample(SCALAR step, SCALAR epsi=SCALAR(0))
	{
		std::cout << "BB: [" << bb_min_.transpose() << " - " << bb_max_.transpose() << "]" << std::endl;
		Epsilon = epsi==SCALAR(0) ? step / SCALAR(20): epsi;
		std::vector<VEC3> samp1;
		samp1.reserve(8192);
		std::vector<VEC3> samp2;
		samp2.reserve(8192);
		std::vector<VEC3> samp3;
		samp3.reserve(8192);

		// Z
		std::thread th1([&]()
		{
		for (SCALAR y = bb_min_[1] + step / 2; y < bb_max_[1] - step / 2; y += step)
			for (SCALAR x = bb_min_[0] + step / 2; x < bb_max_[0] - step / 2; x += step)
				inter_skeleton(VEC3{x, y, bb_min_[2]}, VEC3{0, 0, 1}, samp1, bb_max_[2] - bb_min_[2]);
		});

		// Y
		std::thread th2([&]()
		{
		for (SCALAR x = bb_min_[0] + step / 2; x < bb_max_[0] - step / 2; x += step)
			for (SCALAR z = bb_min_[2] + step / 2; z < bb_max_[2] - step / 2; z += step)
				inter_skeleton(VEC3{x, bb_min_[1], z}, VEC3{0, 1, 0}, samp2, bb_max_[1] - bb_min_[1]);
		});

		// X
		std::thread th3([&]()
		{
		for (SCALAR y = bb_min_[1] + step / 2; y < bb_max_[1] - step / 2; y += step)
			for (SCALAR z = bb_min_[2] + step / 2; z < bb_max_[2] - step / 2; z += step)
				inter_skeleton(VEC3{bb_min_[0], y, z}, VEC3{1, 0, 0}, samp3, bb_max_[0] - bb_min_[0]);
		});
		th1.join();
		th2.join();
		th3.join();

		samples_.clear();
		samples_.insert(samples_.end(),samp1.begin(),samp1.end());
		samples_.insert(samples_.end(),samp2.begin(),samp2.end());
		samples_.insert(samples_.end(),samp3.begin(),samp3.end());

		std::cout << samples_.size() << "points genrated" << std::endl;
	}


	std::vector<VEC3> samples()
	{
		return samples_;
	}

	const VEC3& BBmin() const { return bb_min_;}
	const VEC3& BBmax() const { return bb_max_;}
	VEC3 BBcenter() const { return (bb_min_+bb_max_)/SCALAR(2);}
	VEC3 BBwidth() const { return bb_max_-bb_min_;}

protected:
	std::vector<VEC4> vertices_;  // one by one
	std::vector<VEC4> edges_;	  // by pair
	std::vector<VEC4> triangles_; // by 5

	std::vector<VEC3> samples_;
	VEC3 bb_min_;
	VEC3 bb_max_;

	SCALAR Epsilon;

	inline void update_BB(const VEC4& v)
	{
		// BB
		for (int i = 0; i < 3; ++i)
		{
			SCALAR s = v[i] - v[3];
			if (s < bb_min_[i])
				bb_min_[i] = s;
			s = v[i] + v[3];
			if (s > bb_max_[i])
				bb_max_[i] = s;
		}
	}

	inline void update_BB(const VEC3& v, SCALAR r)
	{
		update_BB(VEC4(v[0], v[1], v[2], r));
	}

	static inline SCALAR evalSphereSDF(const VEC3& P, typename std::vector<VEC4>::const_iterator it)
	{
		return (it->template topRows<3>() - P).norm() - it->w();
	}

	static inline SCALAR evalConeSDF(const VEC3& P, typename std::vector<VEC4>::const_iterator it)
	{
		const VEC3& P_a = it->template topRows<3>();
		SCALAR ra = it->w();
		++it;
		const VEC3& Pb = it->template topRows<3>();
		SCALAR rb = it->w();

		SCALAR rba = rb - ra;
		VEC3 ba = Pb - P_a;
		SCALAR baba = ba.dot(ba);
		VEC3 pa = P - P_a;
		SCALAR papa = pa.dot(pa);
		SCALAR paba = pa.dot(ba) / baba;
		SCALAR x = std::sqrt(papa - paba * paba * baba);
		SCALAR cax = std::max(SCALAR(0), x - ((paba < SCALAR(0.5)) ? ra : rb));
		SCALAR cay = std::abs(paba - SCALAR(0.5)) - SCALAR(0.5);
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
		SCALAR dist = evalPlaneSDF(P, it);
		for (int i = 1; i < 5; ++i)
			dist = std::max(dist, evalPlaneSDF(P, ++it));
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

		compute_CenterDisk(Sph1, Sph2, centers[1], centers[2]);
		VEC3 N1 = (Sph2.template topRows<3>() - Sph1.template topRows<3>()).normalized();
		compute_CenterDisk(Sph2, Sph3, centers[3], centers[4]);
		VEC3 N2 = (Sph3.template topRows<3>() - Sph2.template topRows<3>()).normalized();
		compute_CenterDisk(Sph3, Sph1, centers[5], centers[0]);
		VEC3 N3 = (Sph1.template topRows<3>() - Sph3.template topRows<3>()).normalized();

		std::array<VEC3, 6> I;

		interPPS({N1[0], N1[1], N1[2], -N1.dot(centers[1])}, {-N3[0], -N3[1], -N3[2], N3.dot(centers[0])}, Sph1, I[0],
				 I[3]);
		interPPS({N2[0], N2[1], N2[2], -N2.dot(centers[3])}, {-N1[0], -N1[1], -N1[2], N1.dot(centers[2])}, Sph2, I[1],
				 I[5]);
		interPPS({N3[0], N3[1], N3[2], -N3.dot(centers[5])}, {-N2[0], -N2[1], -N2[2], N2.dot(centers[4])}, Sph3, I[2],
				 I[4]);

		Planes[0].template topRows<3>() = (I[1] - I[0]).cross(I[2] - I[0]).normalized();
		Planes[0].w() = -Planes[0].template topRows<3>().dot(I[0]);

		Planes[1].template topRows<3>() = (I[4] - I[3]).cross(I[5] - I[3]).normalized();
		Planes[1].w() = -Planes[1].template topRows<3>().dot(I[3]);

		Planes[2].template topRows<3>() = (I[5] - I[3]).cross(I[0] - I[3]).normalized();
		Planes[2].w() = -Planes[2].template topRows<3>().dot(I[3]);

		Planes[3].template topRows<3>() = (I[4] - I[5]).cross(I[1] - I[5]).normalized();
		Planes[3].w() = -Planes[3].template topRows<3>().dot(I[5]);

		Planes[4].template topRows<3>() = (I[3] - I[4]).cross(I[2] - I[4]).normalized();
		Planes[4].w() = -Planes[4].template topRows<3>().dot(I[4]);
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_RENDERING_SHAPE3_DRAWER_H_

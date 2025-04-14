#ifndef CGOGN_SIMULATION_XPBD_XPBD_H_
#define CGOGN_SIMULATION_XPBD_XPBD_H_

#include <cgogn/core/types/maps/map_base.h>

#include <algorithm>
#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/volume.h>
#include <cgogn/geometry/types/vector_traits.h>

#define POISSON_RATIO 0.45
#define YOUNG_MODULUS 1e7

#define LAME_MU (YOUNG_MODULUS / (2 * (1 + POISSON_RATIO)))
#define LAME_LAMBDA ((YOUNG_MODULUS * POISSON_RATIO) / ((1 + POISSON_RATIO) * (1 - 2 * POISSON_RATIO)))

#define SHEAR_MODULUS (YOUNG_MODULUS / (2 * (1 + POISSON_RATIO)))
#define BULK_MODULUS (YOUNG_MODULUS / (3 * (1 - 2 * POISSON_RATIO)))

#define NUM_SUBSTEP 5
#define DENSITY 10

#define EPS 1e-12

namespace cgogn
{
namespace simulation
{
template <typename MAP>
class XPBD
{
	using Self = XPBD;
	template <typename T>
	using Attribute = typename mesh_traits<MAP>::template Attribute<T>;
	using Vec3 = geometry::Vec3;
	using Mat3d = geometry::Mat3d;
	using Vertex = typename mesh_traits<MAP>::Vertex;
	using Volume = typename mesh_traits<MAP>::Volume;

public:
	// Stable Values
	std::shared_ptr<Attribute<Vec3>> init_pos_;
	std::shared_ptr<Attribute<Vec3>> init_cm_;
	std::shared_ptr<Attribute<double>> init_volume_;
	std::shared_ptr<Attribute<double>> masse_;
	std::shared_ptr<Attribute<Mat3d>> inv_Q_;
	std::shared_ptr<Attribute<bool>> fixed_vertex;
	std::shared_ptr<Attribute<std::vector<Vertex>>> inc_vertices_;
	// Integration Values

	std::shared_ptr<Attribute<Vec3>> pos_;
	std::shared_ptr<Attribute<Vec3>> pos_prev_;
	std::shared_ptr<Attribute<double>> volume_;
	std::shared_ptr<Attribute<Vec3>> speed_;
	std::shared_ptr<Attribute<Vec3>> f_ext_;
	std::shared_ptr<Attribute<Vec3>> Grad_C_i_;
	std::shared_ptr<Attribute<Vec3>> Grad_C2_i_;

	XPBD()
		: init_pos_(nullptr), init_cm_(nullptr), masse_(nullptr), inc_vertices_(nullptr), pos_(nullptr),
		  pos_prev_(nullptr), speed_(nullptr), f_ext_(nullptr), inv_Q_(nullptr), Grad_C_i_(nullptr), volume_(nullptr)
	{
	}

	void init_solver(MAP& m, std::shared_ptr<Attribute<Vec3>> pos)
	{
		pos_ = pos;
		init_pos_ = add_attribute<Vec3, Vertex>(m, "XPBD_Init_pos");
		init_cm_ = add_attribute<Vec3, Volume>(m, "XPBD_init_cm");

		masse_ = add_attribute<double, Vertex>(m, "XPBD_masse");

		init_volume_ = add_attribute<double, Volume>(m, "XPBD_init_volume");
		inv_Q_ = add_attribute<Mat3d, Volume>(m, "XPBD_inv_Q");
		inc_vertices_ = add_attribute<std::vector<Vertex>, Volume>(m, "XPBD_inc_vertices_vector_");

		pos_prev_ = add_attribute<Vec3, Vertex>(m, "XPBD_pos_prev_");
		speed_ = add_attribute<Vec3, Vertex>(m, "XPBD_speed_");
		volume_ = add_attribute<double, Volume>(m, "XPBD_volume");

		f_ext_ = add_attribute<Vec3, Vertex>(m, "XPBD_f_ext");

		Grad_C_i_ = add_attribute<Vec3, Vertex>(m, "XPBD_Grad_C_i");
		Grad_C2_i_ = add_attribute<Vec3, Vertex>(m, "XPBD_Grad_C2_i");

		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			value<Vec3>(m, init_pos_, v) = value<Vec3>(m, pos_, v);
			value<double>(m, masse_, v) = 0;
			value<Vec3>(m, speed_, v) = Vec3(0, 0, 0);
			value<Vec3>(m, f_ext_, v) = Vec3(0, 0, 0);
			return true;
		});
		
		geometry::compute_volume(m, pos_.get(), init_volume_.get());
		foreach_cell(m, [&](Volume v) -> bool {
			std::vector<Vertex> vertices;
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				vertices.push_back(w);
				return true;
			});
			double masse = value<double>(m, init_volume_, v) * DENSITY / vertices.size();
			for (auto w : vertices)
			{
				value<double>(m, masse_, w) += masse;
			}
			return true;
		});

		parallel_foreach_cell(m, [&](Volume v) -> bool {
			double masse = 0;
			Vec3 cm = Vec3(0, 0, 0);
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				masse += value<double>(m, masse_, w);
				cm += value<double>(m, masse_, w) * value<Vec3>(m, init_pos_, w);
				return true;
			});
			value<Vec3>(m, init_cm_, v) = cm / masse;
			return true;
		});


		parallel_foreach_cell(m, [&](Volume v) -> bool {
			Mat3d Q = Mat3d::Zero();
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				Vec3 init_r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
				double masse = value<double>(m, masse_, w);
				Q += masse * init_r_i * init_r_i.transpose();
				return true;
			});
			value<Mat3d>(m, inv_Q_, v) = Q.inverse();
			return true;
		});
	}

	void constraint_Neo_Hookean_H(MAP& m, Volume v, double h)
	{
		// Compute center of mass
		Vec3 cm = Vec3(0,0,0);
		double masse = 0;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			masse += value<double>(m, masse_, w);
			cm += value<double>(m, masse_, w) * value<Vec3>(m, pos_, w);
			return true;
		});
		cm = cm / masse; 
		
		// Compute P
		Mat3d P = Mat3d::Zero();
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			Vec3 ri = value<Vec3>(m, pos_, w) - cm;
			Vec3 r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
			double masse = value<double>(m, masse_, w);
			P += masse * ri * r_i.transpose();
			return true;
		});

		Mat3d inv_Q = value<Mat3d>(m, inv_Q_, v);
		// compute F = P*Q^-1
        Mat3d F = P * inv_Q;

        //Petit bout de code poour la stabilité numérique ne pas toucher
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (i == j)
				{
					if (fabs(F(i, j) - 1) < EPS)
						F(i, j) = 1;
				}
				else
				{
					if (fabs(F(i, j)) < EPS)
						F(i, j) = 0;
				}
			}
		}

        // Compute C  = det(F) - (1+MU/LAMBDA)

		double C = F.determinant() - (1.+(LAME_MU/LAME_LAMBDA));
        // Compute Volume
        double Ve = 0;
		Ve = cgogn::geometry::volume(m, v , pos_.get());
		
		 
		// Compute alpha_H = 1/(LAMBDA*V)
		double alpha_h = 1./(LAME_LAMBDA*Ve);

		// Compute denum
		double denum = 0;

		Vec3 f23 = F.col(1).cross(F.col(2));
		Vec3 f31 = F.col(2).cross(F.col(0));
		Vec3 f12 = F.col(0).cross(F.col(1));

		Mat3d F_vec = Mat3d();
		F_vec(0,0) = f23[0];
		F_vec(0,1) = f23[1];
		F_vec(0,2) = f23[2];
		F_vec(1,0) = f31[0];
		F_vec(1,1) = f31[1];
		F_vec(1,2) = f31[2];
		F_vec(2,0) = f12[0];
		F_vec(2,1) = f12[1];
		F_vec(2,2) = f12[2];

        foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
            // Compute Grad_C_i
			Vec3 res;
			double masse = value<double>(m,masse_,w);
			Mat3d inv_T_Q = value<Mat3d>(m,inv_Q_,v).transpose();
			Vec3 r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
			res = masse * F_vec * inv_T_Q * r_i;
			value<Vec3>(m,Grad_C_i_,w) = res;
			double norm_res = res.squaredNorm();
			denum += (1./masse) * norm_res ;
            return true;
        });

		// Compute lambda
		denum += alpha_h/(h*h);
        double lambda;
		lambda = (-C) / denum;

        foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
            // Compute and apply delta_x_i
			double masse = value<double>(m,masse_,w);
			value<Vec3>(m,pos_,w) = value<Vec3>(m,pos_,w) + (lambda * (1/masse) * value<Vec3>(m,Grad_C_i_,w));
            return true;
        });
	}
	void constraint_Neo_Hookean_D(MAP& m, Volume v, double h)
	{
		// Compute center of mass

		Vec3 cm = Vec3(0,0,0);
		double masse = 0;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			masse += value<double>(m, masse_, w);
			cm += value<double>(m, masse_, w) * value<Vec3>(m, pos_, w);
			return true;
		});
		cm = cm / masse; 

		// Compute P
		Mat3d P = Mat3d::Zero();
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			Vec3 ri = value<Vec3>(m, pos_, w) - cm;
			Vec3 r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
			double masse = value<double>(m, masse_, w);
			P += masse * ri * r_i.transpose();
			return true;
		});

		Mat3d inv_Q = value<Mat3d>(m, inv_Q_, v);
		// compute F = P*Q^-1
        Mat3d F = P * inv_Q;

		double r = sqrt(F.col(0).squaredNorm() + F.col(1).squaredNorm() + F.col(2).squaredNorm());
        foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
            Vec3 res;
			double masse = value<double>(m,masse_,w);
			Mat3d inv_T_Q = value<Mat3d>(m,inv_Q_,v).transpose();
			Vec3 r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
			res = (masse/r) * F * inv_T_Q * r_i;
			value<Vec3>(m,Grad_C2_i_,w) = res;
            return true;
        });

		//Petit bout de code poour la stabilité numérique ne pas toucher
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (i == j)
				{
					if (fabs(F(i, j) - 1) < EPS)
						F(i, j) = 1;
				}
				else
				{
					if (fabs(F(i, j)) < EPS)
						F(i, j) = 0;
				}
			}
		}

        // Compute C  = sqrt(tr(F^T*F))

		double C = sqrt((F.transpose()*F).trace());

        // Compute Volume

		double Ve = 0;
		Ve = cgogn::geometry::volume(m, v , pos_.get());

        // Compute alpha_D = 1/(MU*V)

		double alpha_D = 1/(LAME_MU*Ve);

		// Compute denum
		double denum = 0;

		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			double masse = value<double>(m,masse_,w);
			double norm_grad = value<Vec3>(m,Grad_C2_i_,w).squaredNorm();
			denum += (1/masse) * norm_grad;
            return true;
        });

		denum += alpha_D/(h*h);

        // Compute lambda

		double lambda;
		lambda = (-C) / denum;

        //Compute and apply delta_x for each vertices

		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			double masse = value<double>(m,masse_,w);
			value<Vec3>(m,pos_,w) = value<Vec3>(m,pos_,w) + (lambda * (1/masse) * value<Vec3>(m,Grad_C2_i_,w));
            return true;
        });

	}
	void constraint_Zero_Energy(MAP& m, Volume v, double)
	{
		// Compute center of mass
		Vec3 cm = Vec3(0,0,0);
		double masse = 0;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			masse += value<double>(m, masse_, w);
			cm += value<double>(m, masse_, w) * value<Vec3>(m, pos_, w);
			return true;
		});
		cm = cm / masse; 

		// Compute P
		Mat3d P = Mat3d::Zero();
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			Vec3 ri = value<Vec3>(m, pos_, w) - cm;
			Vec3 r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
			double masse = value<double>(m, masse_, w);
			P += masse * ri * r_i.transpose();
			return true;
		});

		Mat3d inv_Q = value<Mat3d>(m, inv_Q_, v);
		// compute F = P*Q^-1
        Mat3d F = P * inv_Q;

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (i == j)
				{
					if (fabs(F(i, j) - 1) < EPS)
						F(i, j) = 1;
				}
				else
				{
					if (fabs(F(i, j)) < EPS)
						F(i, j) = 0;
				}
			}
		}
		// apply eq (15)
        foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			Vec3 r_i = value<Vec3>(m, init_pos_, w) - value<Vec3>(m, init_cm_, v);
			value<Vec3>(m,pos_,w) = cm + (F * r_i);
            return true;
        });
		
	}


	void applyDamping(MAP& m, Volume v, double damping_coeff, double time_step)
	{
		Vec3 x_cm = Vec3::Zero();
		Vec3 v_cm = Vec3::Zero();
		Vec3 L = Vec3::Zero();
		Mat3d I = Mat3d::Zero();
		double sm_i = 0;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			double m_i = value<double>(m, masse_, w);
			x_cm += m_i * value<Vec3>(m, pos_.get(), w);
			v_cm += m_i * value<Vec3>(m, speed_, w);
			sm_i += m_i;
			return true;
		});
		x_cm /= sm_i;
		v_cm /= sm_i;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			Vec3 r_i = value<Vec3>(m, pos_.get(), w) - x_cm;
			double m_i = value<double>(m, masse_, w);
			L += r_i.cross(m_i * value<Vec3>(m, speed_, w));
			Mat3d Ri;
			Ri << 0, -r_i(2), r_i(1), r_i(2), 0, -r_i(0), -r_i(1), r_i(0), 0;
			I += Ri * Ri.transpose() * m_i;
			return true;
		});
		Vec3 omega = I.inverse() * L;
		foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
			Vec3 r_i = value<Vec3>(m, pos_.get(), w) - x_cm;
			Vec3 new_v_i = v_cm + omega.cross(r_i);
			value<Vec3>(m, speed_, w) +=
				std::min(damping_coeff * time_step, 1.0) * (new_v_i - value<Vec3>(m, speed_, w));
			return true;
		});
	}

	void solver(MAP& m, double timestep)
	{

		double h = timestep / NUM_SUBSTEP;
		std::vector<Volume> vec_volume;
		std::vector<Vertex> vec_vertices;
		foreach_cell(m, [&](Vertex v) -> bool {
			vec_vertices.push_back(v);
			return true;
		});
		foreach_cell(m, [&](Volume v) -> bool {
			vec_volume.push_back(v);
			value<std::vector<Vertex>>(m, inc_vertices_.get(), v).clear();
			foreach_incident_vertex(m, v, [&](Vertex w) -> bool {
				value<std::vector<Vertex>>(m, inc_vertices_.get(), v).push_back(w);
				return true;
			});
			return true;
        });

		for (int i = 0; i < NUM_SUBSTEP; i++)
		{
			// Initialisation sub step
			for (Vertex v : vec_vertices)
			{
				if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), v))
				{
					continue;
				}
				value<Vec3>(m, pos_prev_, v) = value<Vec3>(m, pos_, v);
				value<Vec3>(m, speed_, v) += h * value<Vec3>(m, f_ext_, v) / value<double>(m, masse_, v);
				value<Vec3>(m, pos_, v) += h * value<Vec3>(m, speed_, v);
			}
			// Constraint solver
			for (Volume v : vec_volume)
			{
				constraint_Neo_Hookean_H(m, v, h);
				constraint_Neo_Hookean_D(m, v, h);
				constraint_Zero_Energy(m, v, h);
			}
			// Speed Solver
			for (Vertex v : vec_vertices)
			{
				if (this->fixed_vertex && value<bool>(m, this->fixed_vertex.get(), v))
				{
					continue;
				}
				Vec3 new_v = (value<Vec3>(m, pos_, v) - value<Vec3>(m, pos_prev_, v)) / h;
				for (int i = 0; i < 3; i++)
					if (fabs(new_v[i]) < EPS)
						new_v[i] = 0;
				// value<Vec3>(m, speed_, v) = (1 - (0.005 * h)) * new_v;
				value<Vec3>(m, speed_, v) = new_v;
			}
			// Damping
			foreach_cell(m, [&](Volume v) -> bool {
				applyDamping(m, v, 0.1, timestep);
				return true;
			});
		}
	}
};
} // namespace simulation
} // namespace cgogn
#endif // CGOGN_SIMULATION_XPBD_XPBD_H_

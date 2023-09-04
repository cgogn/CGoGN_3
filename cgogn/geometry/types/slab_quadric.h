#ifndef CGOGN_GEOMETRY_TYPES_SLAB_QUADRIC_H_
#define CGOGN_GEOMETRY_TYPES_SLAB_QUADRIC_H_
#include <cgogn/geometry/types/vector_traits.h>
#include <cgogn/modeling/algos/decimation/SQEM_helper.h>
namespace cgogn
{

namespace geometry
{
class Slab_Quadric
{

public:
	Slab_Quadric()
	{
	}
 	inline Slab_Quadric(const Vec4& p1, const Vec4& n1, const Vec4& n2,
						bool boundary = false)
	{
		
		if (!boundary)
		{
			_A = (n1 * n1.transpose()) + (n2 * n2.transpose());
			
			_b = -2 * _A * p1 ;
			_c = p1.transpose() * _A * p1 ;
			_add_A.setZero();
			_add_b.setZero();
			_add_c = 0;
		}
		else
		{
			_A.setZero();
			_b.setZero();
			_c = 0;
			_add_A = (n1 * n1.transpose()) + (n2 * n2.transpose());
			_add_b = -2 * _add_A * p1;
			_add_c = p1.transpose() * _add_A * p1;
		}
	}

	inline Slab_Quadric(const Vec4& p1, const Vec4& n1, double stability_ratio)
	{
		_A.setZero();
		_b.setZero();
		_c = 0;
		_add_A = 2* n1 * n1.transpose() * 0.1 * stability_ratio*stability_ratio;
		_add_b = -2 * _add_A * p1;
		_add_c = p1.transpose() * _add_A * p1;
	}
	Slab_Quadric& operator=(const Slab_Quadric& q)
	{
		_A = q._A;
		_b = q._b;
		_c = q._c;
		_add_A = q._add_A;
		_add_b = q._add_b;
		_add_c = q._add_c;
		return *this;
	}

	Slab_Quadric& operator+=(const Slab_Quadric& q)
	{
		_A += q._A;
		_b += q._b;
		_c += q._c;
		_add_A += q._add_A;
		_add_b += q._add_b;
		_add_c += q._add_c;
		return *this;
	}
	
	friend Slab_Quadric operator+(Slab_Quadric lhs, Slab_Quadric& rhs)
	{
		lhs += rhs;
		return lhs;
	}

	inline Scalar eval(const Vec4& center) const
	{
		Scalar cost = center.transpose() * _A* center;
		cost += _b.transpose() * center;
		cost += _c ;
		return cost;
	}

	bool optimized(Vec4& v)
	{
		Mat4 m(_A);
		Mat4 inverse;
		Scalar determinant;
		bool invertible;
		m.computeInverseAndDetWithCheck(inverse, determinant, invertible, 0.01);
		if (invertible)
			v = -0.5* inverse * _b ;
		return invertible;
	}

	void clear(){
		_A.setZero();
		_b.setZero();
		_c = 0;
		_add_A.setZero();
		_add_b.setZero();
		_add_c = 0;
	}


	Mat4 _A;
	Vec4 _b;
	Scalar _c = 0;
	Mat4 _add_A;
	Vec4 _add_b;
	Scalar _add_c;
};
} // namespace geometry
} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_SLAB_QUADRIC_H_
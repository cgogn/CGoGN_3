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

#include <cgogn/geometry/functions/distance.h>

namespace cgogn
{

namespace geometry
{

// http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf

Scalar squared_distance_point_triangle(const Vec3& P, const Vec3& A, const Vec3& B, const Vec3& C)
{
	Vec3 D = A - P;
	Vec3 E0 = B - A;
	Vec3 E1 = C - A;

	Scalar a = E0.squaredNorm();
	Scalar b = E0.dot(E1);
	Scalar c = E1.squaredNorm();
	Scalar d = E0.dot(D);
	Scalar e = E1.dot(D);
	Scalar f = D.squaredNorm();

	Scalar det = fabs(a * c - b * b);
	Scalar s = b * e - c * d;
	Scalar t = b * d - a * e;

	Scalar squared_distance;

	if (s + t <= det)
	{
		if (s < 0.0)
		{
			if (t < 0.0) // region 4
			{
				if (d < 0.0)
				{
					t = 0.0;
					if (-d >= a)
					{
						s = 1.0;
						squared_distance = a + 2.0 * d + f;
					}
					else
					{
						s = -d / a;
						squared_distance = d * s + f;
					}
				}
				else
				{
					s = 0.0;
					if (e >= 0.0)
					{
						t = 0.0;
						squared_distance = f;
					}
					else if (-e >= c)
					{
						t = 1.0;
						squared_distance = c + 2.0 * e + f;
					}
					else
					{
						t = -e / c;
						squared_distance = e * t + f;
					}
				}
			}
			else // region 3
			{
				s = 0.0;
				if (e >= 0.0)
				{
					t = 0.0;
					squared_distance = f;
				}
				else if (-e >= c)
				{
					t = 1.0;
					squared_distance = c + 2.0 * e + f;
				}
				else
				{
					t = -e / c;
					squared_distance = e * t + f;
				}
			}
		}
		else if (t < 0.0) // region 5
		{
			t = 0.0;
			if (d >= 0.0)
			{
				s = 0.0;
				squared_distance = f;
			}
			else if (-d >= a)
			{
				s = 1.0;
				squared_distance = a + 2.0 * d + f;
			}
			else
			{
				s = -d / a;
				squared_distance = d * s + f;
			}
		}
		else // region 0
		{
			// minimum at interior point
			Scalar invDet = 1.0 / det;
			s *= invDet;
			t *= invDet;
			squared_distance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f;
		}
	}
	else
	{
		Scalar tmp0, tmp1, numer, denom;

		if (s < 0.0) // region 2
		{
			tmp0 = b + d;
			tmp1 = c + e;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a - 2.0 * b + c;
				if (numer >= denom)
				{
					s = 1.0;
					t = 0.0;
					squared_distance = a + 2.0 * d + f;
				}
				else
				{
					s = numer / denom;
					t = 1.0 - s;
					squared_distance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f;
				}
			}
			else
			{
				s = 0.0;
				if (tmp1 <= 0.0)
				{
					t = 1.0;
					squared_distance = c + 2.0 * e + f;
				}
				else if (e >= 0.0)
				{
					t = 0.0;
					squared_distance = f;
				}
				else
				{
					t = -e / c;
					squared_distance = e * t + f;
				}
			}
		}
		else if (t < 0.0) // region 6
		{
			tmp0 = b + e;
			tmp1 = a + d;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a - 2.0 * b + c;
				if (numer >= denom)
				{
					t = 1.0;
					s = 0.0;
					squared_distance = c + 2.0 * e + f;
				}
				else
				{
					t = numer / denom;
					s = 1.0 - t;
					squared_distance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f;
				}
			}
			else
			{
				t = 0.0;
				if (tmp1 <= 0.0)
				{
					s = 1.0;
					squared_distance = a + 2.0 * d + f;
				}
				else if (d >= 0.0)
				{
					s = 0.0;
					squared_distance = f;
				}
				else
				{
					s = -d / a;
					squared_distance = d * s + f;
				}
			}
		}
		else // region 1
		{
			numer = c + e - b - d;
			if (numer <= 0.0)
			{
				s = 0.0;
				t = 1.0;
				squared_distance = c + 2.0 * e + f;
			}
			else
			{
				denom = a - 2.0 * b + c;
				if (numer >= denom)
				{
					s = 1.0;
					t = 0.0;
					squared_distance = a + 2.0 * d + f;
				}
				else
				{
					s = numer / denom;
					t = 1.0 - s;
					squared_distance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f;
				}
			}
		}
	}

	//	squared_distance = s * (a*s + b*t + 2.0*d) + t * (b*s + c*t + 2.0*e) + f;

	// account for numerical round-off error
	if (squared_distance < 0.0)
		squared_distance = 0.0;

	return squared_distance;
}

void closest_point_in_triangle(const Vec3& P, const Vec3& A, const Vec3& B, const Vec3& C, Scalar& u, Scalar& v,
							   Scalar& w)
{
	Vec3 D = A - P;
	Vec3 E0 = B - A;
	Vec3 E1 = C - A;

	Scalar a = E0.squaredNorm();
	Scalar b = E0.dot(E1);
	Scalar c = E1.squaredNorm();
	Scalar d = E0.dot(D);
	Scalar e = E1.dot(D);
	// Scalar f = D.squaredNorm() ;

	Scalar det = fabs(a * c - b * b);
	Scalar s = b * e - c * d;
	Scalar t = b * d - a * e;

	if (s + t <= det)
	{
		if (s < 0.0)
		{
			if (t < 0.0) // region 4
			{
				if (d < 0.0)
				{
					t = 0.0;
					if (-d >= a)
						s = 1.0;
					else
						s = -d / a;
				}
				else
				{
					s = 0.0;
					if (e >= 0.0)
						t = 0.0;
					else if (-e >= c)
						t = 1.0;
					else
						t = -e / c;
				}
			}
			else // region 3
			{
				s = 0.0;
				if (e >= 0.0)
					t = 0.0;
				else if (-e >= c)
					t = 1.0;
				else
					t = -e / c;
			}
		}
		else if (t < 0.0) // region 5
		{
			t = 0.0;
			if (d >= 0.0)
				s = 0.0;
			else if (-d >= a)
				s = 1.0;
			else
				s = -d / a;
		}
		else // region 0
		{
			// minimum at interior point
			Scalar invDet = 1.0 / det;
			s *= invDet;
			t *= invDet;
		}
	}
	else
	{
		Scalar tmp0, tmp1, numer, denom;

		if (s < 0.0) // region 2
		{
			tmp0 = b + d;
			tmp1 = c + e;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a - 2.0 * b + c;
				if (numer >= denom)
				{
					s = 1.0;
					t = 0.0;
				}
				else
				{
					s = numer / denom;
					t = 1.0 - s;
				}
			}
			else
			{
				s = 0.0;
				if (tmp1 <= 0.0)
					t = 1.0;
				else if (e >= 0.0)
					t = 0.0;
				else
					t = -e / c;
			}
		}
		else if (t < 0.0) // region 6
		{
			tmp0 = b + e;
			tmp1 = a + d;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a - 2.0 * b + c;
				if (numer >= denom)
				{
					t = 1.0;
					s = 0.0;
				}
				else
				{
					t = numer / denom;
					s = 1.0 - t;
				}
			}
			else
			{
				t = 0.0;
				if (tmp1 <= 0.0)
					s = 1.0;
				else if (d >= 0.0)
					s = 0.0;
				else
					s = -d / a;
			}
		}
		else // region 1
		{
			numer = c + e - b - d;
			if (numer <= 0.0)
			{
				s = 0.0;
				t = 1.0;
			}
			else
			{
				denom = a - 2.0 * b + c;
				if (numer >= denom)
				{
					s = 1.0;
					t = 0.0;
				}
				else
				{
					s = numer / denom;
					t = 1.0 - s;
				}
			}
		}
	}

	//	u = s;
	//	v = t;
	//	w = 1.0 - s - t;

	u = 1.0 - s - t;
	v = s;
	w = t;
}

} // namespace geometry

} // namespace cgogn

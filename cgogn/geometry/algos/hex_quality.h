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

#ifndef CGOGN_GEOMETRY_ALGOS_HEX_QUALITY_H_
#define CGOGN_GEOMETRY_ALGOS_HEX_QUALITY_H_

#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

void compute_hex_frame(const CMap3& m, const CMap3::Attribute<Vec3>* vertex_position,
					   CMap3::Attribute<Mat3>* vertex2_frame, CMap3::Attribute<Mat3>* hex_frame)
{
	parallel_foreach_cell(m, [&](CMap3::Volume v) -> bool {
		Dart d0 = v.dart_;

		Dart D[8];
		D[0] = d0;
		D[1] = phi1(m, d0);
		D[2] = phi1(m, D[1]);
		D[3] = phi1(m, D[2]);
		D[4] = phi<2, 1, 1>(m, d0);
		D[5] = phi<2, 1, 1>(m, D[1]);
		D[6] = phi<2, 1, 1>(m, D[2]);
		D[7] = phi<2, 1, 1>(m, D[3]);

		Vec3 P[8];
		P[0] = value<Vec3>(m, vertex_position, CMap3::Vertex(D[0]));
		P[1] = value<Vec3>(m, vertex_position, CMap3::Vertex(D[1]));
		P[2] = value<Vec3>(m, vertex_position, CMap3::Vertex(D[2]));
		P[3] = value<Vec3>(m, vertex_position, CMap3::Vertex(D[3]));
		P[4] = value<Vec3>(m, vertex_position, CMap3::Vertex(D[4]));
		P[5] = value<Vec3>(m, vertex_position, CMap3::Vertex(D[5]));
		P[6] = value<Vec3>(m, vertex_position, CMap3::Vertex(D[6]));
		P[7] = value<Vec3>(m, vertex_position, CMap3::Vertex(D[7]));

		value<Mat3>(m, hex_frame, v) << ((P[0] + P[1] + P[2] + P[3]) / 4 - (P[4] + P[5] + P[6] + P[7]) / 4),
			((P[0] + P[3] + P[4] + P[7]) / 4 - (P[1] + P[2] + P[5] + P[6]) / 4),
			((P[0] + P[1] + P[4] + P[5]) / 4 - (P[2] + P[3] + P[6] + P[7]) / 4);

		value<Mat3>(m, vertex2_frame, CMap3::Vertex2(D[0])) << (P[1] - P[0]), (P[4] - P[0]), (P[3] - P[0]);
		value<Mat3>(m, vertex2_frame, CMap3::Vertex2(D[1])) << (P[0] - P[1]), (P[2] - P[1]), (P[5] - P[1]);
		value<Mat3>(m, vertex2_frame, CMap3::Vertex2(D[2])) << (P[1] - P[2]), (P[3] - P[2]), (P[6] - P[2]);
		value<Mat3>(m, vertex2_frame, CMap3::Vertex2(D[3])) << (P[0] - P[3]), (P[7] - P[3]), (P[2] - P[3]);
		value<Mat3>(m, vertex2_frame, CMap3::Vertex2(D[4])) << (P[0] - P[4]), (P[5] - P[4]), (P[7] - P[4]);
		value<Mat3>(m, vertex2_frame, CMap3::Vertex2(D[5])) << (P[1] - P[5]), (P[6] - P[5]), (P[4] - P[5]);
		value<Mat3>(m, vertex2_frame, CMap3::Vertex2(D[6])) << (P[2] - P[6]), (P[7] - P[6]), (P[5] - P[6]);
		value<Mat3>(m, vertex2_frame, CMap3::Vertex2(D[7])) << (P[3] - P[7]), (P[4] - P[7]), (P[6] - P[7]);

		return true;
	});
}

void compute_scaled_jacobian(const CMap3& m, const CMap3::Attribute<Mat3>* vertex2_frame,
							 const CMap3::Attribute<Mat3>* hex_frame, CMap3::Attribute<Scalar>* volume_scaled_jacobian)
{
	parallel_foreach_cell(m, [&](CMap3::Volume v) -> bool {
		Mat3 frame_h = value<Mat3>(m, hex_frame, v);
		frame_h.col(0).normalize();
		frame_h.col(1).normalize();
		frame_h.col(2).normalize();
		Scalar jacobian = frame_h.determinant();

		foreach_incident_vertex(m, v, [&](CMap3::Vertex iv) -> bool {
			Mat3 frame = value<Mat3>(m, vertex2_frame, CMap3::Vertex2(iv.dart_));
			frame.col(0).normalize();
			frame.col(1).normalize();
			frame.col(2).normalize();

			Scalar det = frame.determinant();
			jacobian = det < jacobian ? det : jacobian;

			return true;
		});

		value<Scalar>(m, volume_scaled_jacobian, v) = jacobian;

		return true;
	});
}

void compute_jacobian(const CMap3& m, const CMap3::Attribute<Mat3>* vertex2_frame,
					  const CMap3::Attribute<Mat3>* hex_frame, CMap3::Attribute<Scalar>* volume_jacobian)
{
	parallel_foreach_cell(m, [&](CMap3::Volume v) -> bool {
		Mat3 frame_h = value<Mat3>(m, hex_frame, v);
		Scalar jacobian = frame_h.determinant();

		foreach_incident_vertex(m, v, [&](CMap3::Vertex iv) -> bool {
			Mat3 frame = value<Mat3>(m, vertex2_frame, CMap3::Vertex2(iv.dart_));
			frame.col(0).normalize();
			frame.col(1).normalize();
			frame.col(2).normalize();

			Scalar det = frame.determinant();
			jacobian = det < jacobian ? det : jacobian;

			return true;
		});

		value<Scalar>(m, volume_jacobian, v) = jacobian;

		return true;
	});
}

// TODO: go to function
Scalar frame_frobenius(Mat3 frame)
{
	Scalar det = frame.determinant();
	if (det <= std::numeric_limits<Scalar>::min())
		return std::numeric_limits<Scalar>::max();

	Vec3 c0 = frame.col(0);
	Vec3 c1 = frame.col(1);
	Vec3 c2 = frame.col(2);

	Scalar t1 = c0.dot(c0) + c1.dot(c1) + c2.dot(c2);
	Scalar t2 = (c0.cross(c1)).dot(c0.cross(c1)) + (c1.cross(c2)).dot(c1.cross(c2)) + (c2.cross(c0)).dot(c2.cross(c0));

	return sqrt(t1 * t2) / (3 * det);
}

void compute_maximum_aspect_frobenius(const CMap3& m, const CMap3::Attribute<Mat3>* vertex2_frame,
									  CMap3::Attribute<Scalar>* volume_max_froebnius)
{
	foreach_cell(m, [&](CMap3::Volume v) -> bool {
		Scalar frobenius = std::numeric_limits<Scalar>::min(); // TODO: check?

		foreach_incident_vertex(m, v, [&](CMap3::Vertex iv) -> bool {
			Mat3 frame = value<Mat3>(m, vertex2_frame, CMap3::Vertex2(iv.dart_));
			frame.col(0).normalize();
			frame.col(1).normalize();
			frame.col(2).normalize();

			Scalar n = frame_frobenius(frame);
			frobenius = n > frobenius ? n : frobenius;

			return true;
		});

		value<Scalar>(m, volume_max_froebnius, v) = frobenius;

		return true;
	});
}

void compute_mean_aspect_frobenius(const CMap3& m, const CMap3::Attribute<Mat3>* vertex2_frame,
								   CMap3::Attribute<Scalar>* volume_mean_froebnius)
{
	foreach_cell(m, [&](CMap3::Volume v) -> bool {
		Scalar frobenius = 0;

		foreach_incident_vertex(m, v, [&](CMap3::Vertex iv) -> bool {
			Mat3 frame = value<Mat3>(m, vertex2_frame, CMap3::Vertex2(iv.dart_));
			frame.col(0).normalize();
			frame.col(1).normalize();
			frame.col(2).normalize();

			Scalar temp = frame_frobenius(frame);
			frobenius += temp;

			return true;
		});

		value<Scalar>(m, volume_mean_froebnius, v) = frobenius / 8.0;

		return true;
	});
}

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_ALGOS_HEX_QUALITY_H_

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

#ifndef CGOGN_GEOMETRY_TYPES_GRID_H_
#define CGOGN_GEOMETRY_TYPES_GRID_H_

#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/global.h>

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace geometry
{

template <uint32 I, uint32 J, uint32 K, typename MESH>
class Grid
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

public:
	Grid(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
		: mesh_(m), vertex_position_(vertex_position), nb_registered_faces_(0)
	{
		for (uint32 i = 0; i < 3; ++i)
		{
			bb_min_[i] = std::numeric_limits<float64>::max();
			bb_max_[i] = std::numeric_limits<float64>::lowest();
		}
		for (const Vec3& p : *vertex_position_)
		{
			for (uint32 i = 0; i < 3; ++i)
			{
				if (p[i] < bb_min_[i])
					bb_min_[i] = p[i];
				if (p[i] > bb_max_[i])
					bb_max_[i] = p[i];
			}
		}
		Vec3 center = (bb_min_ + bb_max_) / Scalar(2);
		bb_min_ = ((bb_min_ - center) * 1.1) + center;
		bb_max_ = ((bb_max_ - center) * 1.1) + center;
		bb_size_ = bb_max_ - bb_min_;

		foreach_cell(mesh_, [&](Face f) {
			++nb_registered_faces_;
			Vec3 face_bb_min, face_bb_max;
			for (uint32 i = 0; i < 3; ++i)
			{
				face_bb_min[i] = std::numeric_limits<float64>::max();
				face_bb_max[i] = std::numeric_limits<float64>::lowest();
			}
			foreach_incident_vertex(mesh_, f, [&](Vertex v) {
				const Vec3& p = value<Vec3>(mesh_, vertex_position_, v);
				for (uint32 i = 0; i < 3; ++i)
				{
					if (p[i] < face_bb_min[i])
						face_bb_min[i] = p[i];
					if (p[i] > face_bb_max[i])
						face_bb_max[i] = p[i];
				}
				return true;
			});
			Vec3i fmin = point_cell_coord(face_bb_min);
			Vec3i fmax = point_cell_coord(face_bb_max);
			for (uint32 i = fmin[0]; i <= fmax[0]; ++i)
				for (uint32 j = fmin[1]; j <= fmax[1]; ++j)
					for (uint32 k = fmin[2]; k <= fmax[2]; ++k)
						cells_[i][j][k].push_back(f);
			return true;
		});
	}

	Vec3i point_cell_coord(const Vec3& p) const
	{
		Vec3 coord = ((p - bb_min_).array() / bb_size_.array());
		return {int(coord[0] * I), int(coord[1] * J), int(coord[2] * K)};
	}

	template <typename FUNC>
	void foreach_face_at(const Vec3& p, const FUNC& func) const
	{
		static_assert(is_func_parameter_same<FUNC, Face>::value, "Wrong function parameter type");
		Vec3i coord = point_cell_coord(p);
		for (Face f : cells_[coord[0]][coord[1]][coord[2]])
			func(f);
	}

	template <typename FUNC>
	void foreach_face_around(const Vec3& p, const FUNC& func) const
	{
		static_assert(is_func_parameter_same<FUNC, Face>::value, "Wrong function parameter type");
		if (nb_registered_faces_ == 0)
			return;
		Vec3i coord = point_cell_coord(p);
		bool all_empty = true;
		std::array<int, 3> start{coord[0] > 0 ? coord[0] - 1 : 0, coord[1] > 0 ? coord[1] - 1 : 0,
								 coord[2] > 0 ? coord[2] - 1 : 0};
		std::array<int, 3> end{coord[0] < I - 1 ? coord[0] + 1 : int(I - 1),
							   coord[1] < J - 1 ? coord[1] + 1 : int(J - 1),
							   coord[2] < K - 1 ? coord[2] + 1 : int(K - 1)};
		for (int i = start[0]; i <= end[0]; ++i)
			for (int j = start[1]; j <= end[1]; ++j)
				for (int k = start[2]; k <= end[2]; ++k)
				{
					if (cells_[i][j][k].size() > 0)
						all_empty = false;
					for (Face f : cells_[i][j][k])
						func(f);
				}

		while (all_empty)
		{
			start = {start[0] > 0 ? start[0] - 1 : 0, start[1] > 0 ? start[1] - 1 : 0, start[2] > 0 ? start[2] - 1 : 0};
			end = {end[0] < I - 1 ? end[0] + 1 : int(I - 1), end[1] < J - 1 ? end[1] + 1 : int(J - 1),
				   end[2] < K - 1 ? end[2] + 1 : int(K - 1)};
			for (int i = start[0]; i <= end[0]; ++i)
				for (int j = start[1]; j <= end[1]; ++j)
					for (int k = start[2]; k <= end[2]; ++k)
					{
						if (i == start[0] || i == end[0] || j == start[1] || j == end[1] || k == start[2] ||
							k == end[2])
						{
							if (cells_[i][j][k].size() > 0)
								all_empty = false;
							for (Face f : cells_[i][j][k])
								func(f);
						}
					}
		}
	}

private:
	const MESH& mesh_;
	std::shared_ptr<Attribute<Vec3>> vertex_position_;
	Vec3 bb_min_, bb_max_, bb_size_;
	std::vector<Face> cells_[I][J][K];
	uint32 nb_registered_faces_;
};

} // namespace geometry

} // namespace cgogn

#endif // CGOGN_GEOMETRY_TYPES_GRID_H_

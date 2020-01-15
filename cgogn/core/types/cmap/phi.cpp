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

#include <cgogn/core/types/cmap/cph3.h>
#include <cgogn/core/types/cmap/phi.h>

namespace cgogn
{

Dart phi2bis(const CPH3& m, Dart d)
{
	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);

	uint32 face_id = m.face_id(d);
	Dart it = d;

	it = phi2(map, it);

	if (m.face_id(it) == face_id)
		return it;
	else
	{
		do
		{
			it = phi2(map, phi3(map, it));
		} while (m.face_id(it) != face_id);

		return it;
	}
}

Dart phi1(const CPH3& m, Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);

	if (m.current_level_ == m.maximum_level_)
		return phi1(map, d);

	bool finished = false;
	uint32 edge_id = m.edge_id(d);
	Dart it = d;
	do
	{
		it = phi1(map, it);
		if (m.dart_level(it) <= m.current_level_)
			finished = true;
		else
			while (m.edge_id(it) != edge_id)
				it = phi1(map, phi2bis(m, it));
	} while (!finished);
	return it;
}

Dart phi_1(const CPH3& m, Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);

	bool finished = false;
	Dart it = phi_1(map, d);
	uint32 edge_id = m.edge_id(it);

	do
	{
		if (m.dart_level(it) <= m.current_level_)
			finished = true;
		else
		{
			it = phi_1(map, it);
			while (m.edge_id(it) != edge_id)
				it = phi_1(map, phi2bis(m, it));
		}
	} while (!finished);
	return it;
}

Dart phi2(const CPH3& m, Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);
	return phi2(map, phi_1(map, phi1(m, d)));
}

Dart phi3(const CPH3& m, Dart d)
{
	cgogn_message_assert(m.dart_level(d) <= m.current_level_, "Access to a dart introduced after current level");

	const CPH3::CMAP& map = static_cast<const CPH3::CMAP&>(m);
	if (phi3(map, d) == d)
		return d;
	return phi3(map, phi_1(map, phi1(m, d)));
}

} // namespace cgogn

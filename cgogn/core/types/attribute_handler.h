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

#ifndef CGOGN_CORE_ATTRIBUTE_HANDLER_H_
#define CGOGN_CORE_ATTRIBUTE_HANDLER_H_

#include <memory>

namespace cgogn
{

template <typename CELL, typename T, template <typename> class ATTRIBUTE, typename MESH>
class AttributeHandler
{
	MESH* m_;
	std::shared_ptr<ATTRIBUTE<T>> attribute_;

public:
	AttributeHandler(MESH* m, const std::shared_ptr<ATTRIBUTE<T>>& attribute) : m_(m), attribute_(attribute)
	{
	}

	CGOGN_NOT_COPYABLE_NOR_MOVABLE(AttributeHandler);

	inline T& operator[](CELL c)
	{
		return value<T>(*m_, attribute_, c);
	}

	inline const T& operator[](CELL c) const
	{
		return value<T>(*m_, attribute_, c);
	}
};

template <typename CELL, typename T, template <typename> class ATTRIBUTE, typename MESH>
inline AttributeHandler<CELL, T, ATTRIBUTE, MESH> attribute_handler(MESH* m, const std::shared_ptr<ATTRIBUTE<T>>& att)
{
	return AttributeHandler<CELL, T, ATTRIBUTE, MESH>(m, att);
}

} // namespace cgogn

#endif // CGOGN_CORE_ATTRIBUTE_HANDLER_H_

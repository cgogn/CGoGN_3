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

#ifndef CGOGN_SIMULATION_SHALLOW_WATER_RIEMANN_SOLVER_H_
#define CGOGN_SIMULATION_SHALLOW_WATER_RIEMANN_SOLVER_H_

#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace simulation
{

namespace shallow_water
{

using Scalar = geometry::Scalar;

enum BoundaryCondition
{
	BC_F = 0,
	BC_C,
	BC_H,
	BC_Q,
	BC_Z,
	BC_S
};

inline std::string bc_name(BoundaryCondition bc)
{
	switch (bc)
	{
	case BC_F:
		return "Free Outflow";
		break;
	case BC_C:
		return "Critical";
		break;
	case BC_H:
		return "Prescribed H";
		break;
	case BC_Q:
		return "Prescribed Q";
		break;
	case BC_Z:
		return "Prescribed Z";
		break;
	case BC_S:
		return "Weir";
		break;
	default:
		return "";
		break;
	}
}

struct Str_Riemann_Flux
{
	Scalar F1; /**< Flux de masse à travers l'interface **/
	Scalar F2; /**< Flux de quantité de mouvement à travers l'interface dans la direction normale à l'interface **/
	Scalar
		F3; /**< Flux de quantité de mouvement à travers l'interface dans la direction longitudinale à l'interface **/
	Scalar
		s2L; /**< Quantité de mouvement associée well-balancing du terme source pour la maille gauche de l'interface **/
	Scalar
		s2R; /**< Quantité de mouvement associée well-balancing du terme source pour la maille droite de l'interface **/
};

Str_Riemann_Flux Solv_HLLC(Scalar g, Scalar hmin, Scalar smalll, Scalar zbL, Scalar zbR, Scalar PhiL, Scalar PhiR,
						   Scalar hL, Scalar qL, Scalar rL, Scalar hR, Scalar qR, Scalar rR);

Str_Riemann_Flux border_condition(BoundaryCondition typBC, Scalar valBC, Scalar NormX, Scalar NormY, Scalar q, Scalar r,
								  Scalar z, Scalar zb, Scalar g, Scalar hmin, Scalar smalll);

} // namespace shallow_water

} // namespace simulation

} // namespace cgogn

#endif // CGOGN_SIMULATION_SHALLOW_WATER_RIEMANN_SOLVER_H_

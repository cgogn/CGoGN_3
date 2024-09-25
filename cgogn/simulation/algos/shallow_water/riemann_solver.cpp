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

#include <cgogn/simulation/algos/shallow_water/riemann_solver.h>

namespace cgogn
{

namespace simulation
{

namespace shallow_water
{

Str_Riemann_Flux Solv_HLLC(Scalar g, Scalar hmin, Scalar smalll, Scalar zbL, Scalar zbR, Scalar PhiL, Scalar PhiR,
						   Scalar hL, Scalar qL, Scalar rL, Scalar hR, Scalar qR, Scalar rR)
{
	Str_Riemann_Flux Riemann_flux;

	Scalar F1 = 0;
	Scalar F2 = 0;
	Scalar F3 = 0;
	Scalar s2L = 0;
	Scalar s2R = 0;

	Scalar zL = zbL + hL;
	Scalar zR = zbR + hR;

	if (((hL > hmin) && (hR > hmin)) || ((hL < hmin) && (zR >= zbL + hmin) && (hR > hmin)) ||
		((hR < hmin) && (zL >= zbR + hmin) && (hL > hmin)))
	{
		//---possible exchange--------
		// There is water in both cells or one of the cells
		// can fill the other one
		//-----wave speed--------
		Scalar L1L = qL / std::max(hL, smalll) - sqrt(g * std::max(hL, smalll));
		Scalar L3L = qL / std::max(hL, smalll) + sqrt(g * std::max(hL, smalll));
		Scalar L1R = qR / std::max(hR, smalll) - sqrt(g * std::max(hR, smalll));
		Scalar L3R = qR / std::max(hR, smalll) + sqrt(g * std::max(hR, smalll));
		Scalar L1LR = std::min({L1L, L1R, 0e0});
		Scalar L3LR = std::max({L3L, L3R, 0e0});
		//========================
		Scalar PhiLR = std::min(PhiL, PhiR);
		//------compute F1--------
		F1 = L3LR * qL - L1LR * qR + L1LR * L3LR * (zR - zL);
		F1 = F1 * PhiLR / std::max(L3LR - L1LR, smalll);
		//========================
		//-----compute F2---------
		Scalar F2L = (qL * qL) / std::max(hL, smalll) + 5e-1 * g * hL * hL;
		Scalar F2R = (qR * qR) / std::max(hR, smalll) + 5e-1 * g * hR * hR;
		F2 = (L3LR * PhiL * F2L - L1LR * PhiR * F2R + L1LR * L3LR * (PhiR * qR - PhiL * qL)) /
			 std::max(L3LR - L1LR, smalll);
		//==========================
		//-----Compute S2L and S2R---
		Scalar Fact = 0.5 * PhiLR * (hL + hR);
		s2L = 0.5 * (PhiL * hL * hL - PhiR * hR * hR) - Fact * (zL - zR);
		s2L = g * L1LR * s2L / std::max(L3LR - L1LR, smalll);
		s2R = 0.5 * (PhiR * hR * hR - PhiL * hL * hL) - Fact * (zR - zL);
		s2R = g * L3LR * s2R / std::max(L3LR - L1LR, smalll);
		//============================
		//------Compute F3------------
		if (F1 > 0)
			F3 = F1 * rL / std::max(hL, smalll);
		else
			F3 = F1 * rR / std::max(hR, smalll);
	}
	//===================================
	else if ((hL < hmin) && (zR < zbL) && (hR > hmin))
	//------impossible exchange-Cell L empty------
	// The cell L is empty and the water level in the cell R
	// is below zbL-filling is impossible
	{
		F1 = 0e0;
		F2 = 0e0;
		s2L = 0e0;
		s2R = PhiR * 5e-1 * g * hR * hR;
		F3 = 0e0;
	}
	//===============================================
	else if ((hR < hmin) && (zL < zbR) && (hL > hmin))
	//------impossible exchange-Cell R empty------
	// The cell R is empty and the water level in the cell L
	// is below zbR-filling is impossible
	{
		F1 = 0e0;
		F2 = 0e0;
		s2L = -PhiL * 0.5 * g * hL * hL;
		s2R = 0e0;
		F3 = 0e0;
	}
	//===============================================
	else
	// Both cells below hmin:exchange is impossible
	{
		F1 = 0e0;
		F2 = 0e0;
		F3 = 0e0;
		s2L = 0e0;
		s2R = 0e0;
	}

	Riemann_flux.F1 = F1;
	Riemann_flux.F2 = F2;
	Riemann_flux.F3 = F3;
	Riemann_flux.s2L = s2L;
	Riemann_flux.s2R = s2R;

	return Riemann_flux;
}

Str_Riemann_Flux border_condition(BoundaryCondition typBC, Scalar valBC, Scalar NormX, Scalar NormY, Scalar q, Scalar r,
								  Scalar z, Scalar zb, Scalar g, Scalar hmin, Scalar smalll)
{
	Str_Riemann_Flux Flux;

	//-----------initialization------------
	//   h1,q1,r1;h,q,r within the domain
	Scalar q1 = q * NormX + r * NormY;
	Scalar r1 = -q * NormY + r * NormX;
	Scalar h1 = z - zb;

	q1 = -q1;
	r1 = -r1;

	if (h1 < hmin)
	{
		h1 = 0.;
		q1 = 0.;
		r1 = 0.;
	}

	// Characteristic variables
	Scalar c1 = sqrt(g * h1);
	Scalar u1 = q1 / std::max(h1, smalll);
	Scalar v1 = r1 / std::max(h1, smalll);
	Scalar L1 = std::max(u1 + c1, 0.);
	//===================================================================

	Scalar F1 = 0.;
	Scalar F2 = 0.;
	Scalar F3 = 0.;
	//-----Boundary conditions-------------------
	//-----Free Outflow-------
	if (typBC == BC_F)
	{
		//----message------
		F1 = q1;
		F2 = q1 * u1 + 0.5 * g * h1 * h1;
	}
	//=========================
	//-------Critical Section----
	else if (typBC == BC_C)
	{
		//-----message------
		Scalar c = (u1 - 2 * c1) / (valBC - 2);
		c = std::max(c, 0.);
		Scalar u = -valBC * c;
		Scalar h = c * c / g;
		F1 = h * u;
		F2 = h * u * u + 0.5 * g * h * h;
	}
	//============================
	//-------Prescribed h---------
	else if (typBC == BC_H)
	{
		//------message----------
		Scalar h = std::max(valBC, 0.);
		Scalar u = 0;

		if (L1 < 0)
		{
			/* torrentiel sortant*/
			h = h1;
			u = u1;
		}
		else
		{
			Scalar cmin = std::max({sqrt(g * h), (2 * c1 - u1) / 3.0, 0.0});
			h = std::max(h, (cmin * cmin) / g);
			Scalar c = sqrt(g * h);
			u = u1 + 2 * (c - c1);
		}
		F1 = h * u;
		F2 = h * u * u + 0.5 * g * h * h;
	}
	//==============================
	//-------Prescribed z-----------
	else if (typBC == BC_Z)
	{
		//------message-----
		Scalar h = std::max(valBC - zb, 0.);
		Scalar c = sqrt(g * h);
		Scalar u = u1 + 2 * (c - c1);
		F1 = h * u;
		F2 = h * u * u + 0.5 * g * h * h;

		/** @todo Utilité ??? **/
		h = std::max(valBC - zb, 0.); // why is part need
		if (L1 < 0)
		{
			/* torrentiel sortant*/
			h = h1;
			u = u1;
		}
		else
		{
			Scalar cmin = std::max({sqrt(g * h), (2 * c1 - u1) / 3, 0.});
			h = std::max(h, (cmin * cmin) / g);
			c = sqrt(g * h);
			u = u1 + 2 * (c - c1);
		}
		F1 = h * u;
		F2 = h * u * u + 0.5 * g * h * h;
	}
	//===============================
	//--------Prescribed q-----------
	else if (typBC == BC_Q)
	{
		//-----message-------
		F1 = valBC;
		Scalar hc = pow(((F1 * F1) / g), 1 / 3);
		if (hc >= h1)
		{
			F2 = (q1 * q1) / std::max(hc, smalll) + 0.5 * g * h1 * h1 + (F1 - q1) * L1;
		}
		else
		{
			F2 = (q1 * q1) / std::max(h1, smalll) + 0.5 * g * h1 * h1 + (F1 - q1) * L1;
		}
	}
	//=================================
	//---------Weir--------------------
	else if (typBC == BC_S)
	{
		/**
		 ** @todo Implémenter les BC de type 's' en renseignant la cote de la pelle et non la hauteur (permet de gérer
		 *les cas avec plusieurs mailles de cote du fond diférentes attenantes au même seuil)
		 **/
		//-----message-------
		if (h1 < valBC)
		{
			// No z:weir elevation not reached
			F1 = 0;
			F2 = (q1 * q1) / std::max(h1, smalll) + 0.5 * g * h1 * h1;
		}
		else
		{
			// Weir overtoped
			F1 = -0.42 * sqrt(2 * g) * pow((h1 - valBC), 3 / 2);
			F2 = (q1 * q1) / std::max(h1, smalll) + 0.5 * g * h1 * h1;
		}
	}
	else
	{
		// std::cout << "pbl bc" << std::endl;
		// std::cout << typBC << std::endl;
	}

	F3 = (F1 - fabs(F1)) * v1 / 2;

	//----output-----
	F1 = -F1;

	//--return F1,F2,F3
	Flux.F1 = F1;
	Flux.F2 = F2;
	Flux.F3 = F3;
	Flux.s2L = 0.;
	Flux.s2R = 0.;

	return Flux;
}

} // namespace shallow_water

} // namespace simulation

} // namespace cgogn

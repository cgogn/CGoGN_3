/*******************************************************************************
* CGoGN                                                                        *
* Copyright (C) 2019, IGG Group, ICube, University of Strasbourg, France       *
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

#include <cgogn/ui/modules/shallow_water/solver.h>

namespace cgogn
{

namespace ui
{

Str_Riemann_Flux Solv_HLLC(
    Scalar g, Scalar hmin, Scalar smalll,
    Scalar zbL, Scalar zbR,
    Scalar PhiL, Scalar PhiR,
    Scalar hL, Scalar qL, Scalar rL, Scalar hR, Scalar qR, Scalar rR
)
{
    Str_Riemann_Flux Riemann_flux;

	Scalar F1 = 0;
	Scalar F2 = 0;
	Scalar F3 = 0;
	Scalar s2L = 0;
	Scalar s2R = 0;

	Scalar zL = zbL + hL;
	Scalar zR = zbR + hR;

	if (((hL > hmin) && (hR > hmin)) ||
		((hL < hmin) && (zR >= zbL + hmin) && (hR > hmin)) ||
		((hR < hmin) && (zL >= zbR + hmin) && (hL > hmin)))
	{
		//---possible exchange--------
		//There is water in both cells or one of the cells
		//can fill the other one
		//-----wave speed--------
		Scalar L1L = qL / std::max(hL, smalll) - sqrt(g * std::max(hL, smalll));
		Scalar L3L = qL / std::max(hL, smalll) + sqrt(g * std::max(hL, smalll));
		Scalar L1R = qR / std::max(hR, smalll) - sqrt(g * std::max(hR, smalll));
		Scalar L3R = qR / std::max(hR, smalll) + sqrt(g * std::max(hR, smalll));
		Scalar L1LR = std::min({ L1L, L1R, 0e0 });
		Scalar L3LR = std::max({ L3L, L3R, 0e0 });
		//========================
		Scalar PhiLR = std::min(PhiL, PhiR);
		//------compute F1--------
		F1 = L3LR * qL - L1LR * qR + L1LR * L3LR * (zR - zL);
		F1 = F1 * PhiLR / std::max(L3LR - L1LR, smalll);
		//========================
		//-----compute F2---------
		Scalar F2L = (qL * qL) / std::max(hL, smalll) + 5e-1 * g * hL * hL;
		Scalar F2R = (qR * qR) / std::max(hR, smalll) + 5e-1 * g * hR * hR;
		F2 = (L3LR * PhiL * F2L - L1LR * PhiR * F2R + L1LR * L3LR * (PhiR * qR - PhiL * qL))
				/ std::max(L3LR-L1LR, smalll);
		//==========================
		//-----Compute S2L and S2R---
		Scalar Fact = 0.5 * PhiLR * (hL + hR);
		s2L = 0.5 * (PhiL * hL * hL - PhiR * hR * hR)
				- Fact * (zL - zR);
		s2L = g * L1LR * s2L / std::max(L3LR - L1LR, smalll);
		s2R = 0.5 * (PhiR * hR * hR - PhiL * hL * hL)
				- Fact * (zR - zL);
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
		//The cell L is empty and the water level in the cell R
		//is below zbL-filling is impossible
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
		//The cell R is empty and the water level in the cell L
		//is below zbR-filling is impossible
	{
		F1 = 0e0;
		F2 = 0e0;
		s2L = -PhiL * 0.5 * g * hL * hL;
		s2R = 0e0;
		F3 = 0e0;
	}
	//===============================================
	else
		//Both cells below hmin:exchange is impossible
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

} // namespace ui

} // namespace cgogn

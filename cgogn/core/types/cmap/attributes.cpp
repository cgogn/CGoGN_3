/*******************************************************************************
* CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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

#include <cgogn/core/types/cmap/attributes.h>

namespace cgogn
{

AttributeGen::AttributeGen(AttributeContainer* container, bool is_mark, const std::string& name) :
	name_(name),
	container_(container),
	is_mark_(is_mark)
{}

AttributeGen::~AttributeGen()
{
	if (is_mark_)
	{
		auto iter = std::find(container_->mark_attributes_.begin(), container_->mark_attributes_.end(), this);
		if (iter != container_->mark_attributes_.end())
		{
			*iter = container_->mark_attributes_.back();
			container_->mark_attributes_.pop_back();
		}
	}
	else
	{
		auto iter = std::find(container_->attributes_.begin(), container_->attributes_.end(), this);
		if (iter != container_->attributes_.end())
		{
			*iter = container_->attributes_.back();
			container_->attributes_.pop_back();
		}
	}
}

} // namespace cgogn

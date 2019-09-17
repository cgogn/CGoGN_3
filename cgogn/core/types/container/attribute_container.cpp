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

#include <cgogn/core/types/container/attribute_container.h>

#include <cgogn/core/utils/thread_pool.h>
#include <cgogn/core/utils/assert.h>

namespace cgogn
{

/////////////////////////
// AttributeGenT class //
/////////////////////////

AttributeGenT::AttributeGenT(AttributeContainerGen* container, bool is_mark, const std::string& name) :
	container_(container),
	is_mark_(is_mark),
	name_(name)
{}

AttributeGenT::~AttributeGenT()
{
	container_->delete_attribute(this);
}

uint32 AttributeGenT::maximum_index() const
{
	return container_->maximum_index();
}

/////////////////////////////////
// AttributeContainerGen class //
/////////////////////////////////

AttributeContainerGen::AttributeContainerGen() :
	nb_elements_(0),
	maximum_index_(0)
{
	attributes_.reserve(32);
	attributes_shared_ptr_.reserve(32);

	uint32 max_nb_threads = thread_pool()->max_nb_threads();
	mark_attributes_.resize(max_nb_threads);
	available_mark_attributes_.resize(max_nb_threads);
	for (uint32 i = 0; i < max_nb_threads; ++i)
	{
		mark_attributes_[i].reserve(32);
		available_mark_attributes_[i].reserve(32);
	}
	
	available_indices_.reserve(1024);
}

AttributeContainerGen::~AttributeContainerGen()
{}

uint32 AttributeContainerGen::new_index()
{
	uint32 index;
	if (available_indices_.size() > 0)
	{
		index = available_indices_.back();
		available_indices_.pop_back();
	}
	else
		index = maximum_index_++;

	for (AttributeGenT* ag : attributes_)
		ag->manage_index(index);
	
	std::lock_guard<std::mutex> lock(mark_attributes_mutex_);
	for (uint32 i = 0, nb = mark_attributes_.size(); i < nb; ++i)
	{
		for (AttributeGenT* ag : mark_attributes_[i])
			ag->manage_index(index);
	}

	init_ref_counter(index);
	init_mark_attributes(index);
	++nb_elements_;

	return index;
}

void AttributeContainerGen::release_index(uint32 index)
{
	cgogn_message_assert(nb_refs(index) > 0, "Trying to release an unused index");
	available_indices_.push_back(index);
	reset_ref_counter(index);
	--nb_elements_;
}

void AttributeContainerGen::remove_attribute(const std::shared_ptr<AttributeGenT>& attribute)
{
	auto it = std::find(attributes_shared_ptr_.begin(), attributes_shared_ptr_.end(), attribute);
	if (it != attributes_shared_ptr_.end())
	{
		*it = attributes_shared_ptr_.back();
		attributes_shared_ptr_.pop_back();
	}
}

void AttributeContainerGen::remove_attribute(AttributeGenT* attribute)
{
	auto it = std::find_if(
		attributes_shared_ptr_.begin(),
		attributes_shared_ptr_.end(),
		[&] (const auto& att) { return att.get() == attribute; }
	);
	if (it != attributes_shared_ptr_.end())
	{
		*it = attributes_shared_ptr_.back();
		attributes_shared_ptr_.pop_back();
	}
}

void AttributeContainerGen::delete_attribute(AttributeGenT* attribute)
{
	if (attribute->is_mark())
		cgogn_assert_not_reached("Deleting a mark attribute..");
	else
	{
		auto iter = std::find(attributes_.begin(), attributes_.end(), attribute);
		if (iter != attributes_.end())
		{
			*iter = attributes_.back();
			attributes_.pop_back();
		}
	}
}

} // namespace cgogn

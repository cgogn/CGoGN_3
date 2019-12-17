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

#ifndef CGOGN_CORE_CONTAINER_ATTRIBUTE_CONTAINER_H_
#define CGOGN_CORE_CONTAINER_ATTRIBUTE_CONTAINER_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/utils/thread.h>

#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace cgogn
{

class AttributeContainerGen;
template <template <typename> class AttributeT>
class AttributeContainerT;

/////////////////////////
// AttributeGenT class //
/////////////////////////

class CGOGN_CORE_EXPORT AttributeGenT
{
public:
	AttributeGenT(AttributeContainerGen* container, const std::string& name);
	virtual ~AttributeGenT();

	inline const std::string& name() const
	{
		return name_;
	}

	uint32 maximum_index() const;

protected:
	AttributeContainerGen* container_;
	std::string name_;

private:
	friend AttributeContainerGen;
	template <template <typename> class AttributeT>
	friend class AttributeContainerT;

	virtual void manage_index(uint32 index) = 0;
};

/////////////////////////////////
// AttributeContainerGen class //
/////////////////////////////////

class CGOGN_CORE_EXPORT AttributeContainerGen
{
public:
	AttributeContainerGen();
	virtual ~AttributeContainerGen();

	inline uint32 nb_elements() const
	{
		return nb_elements_;
	}
	inline uint32 maximum_index() const
	{
		return maximum_index_;
	}

	uint32 new_index();
	void release_index(uint32 index);

	void remove_attribute(const std::shared_ptr<AttributeGenT>& attribute);
	void remove_attribute(AttributeGenT* attribute);

	using const_iterator = std::vector<std::shared_ptr<AttributeGenT>>::const_iterator;
	inline const_iterator begin() const
	{
		return attributes_shared_ptr_.begin();
	}
	inline const_iterator end() const
	{
		return attributes_shared_ptr_.end();
	}

	inline uint32 first_index() const
	{
		uint32 index = 0u;
		while (index < maximum_index_ && nb_refs(index) == 0)
			++index;
		return index;
	}

	inline uint32 last_index() const
	{
		return maximum_index_;
	}

	inline uint32 next_index(uint32 index) const
	{
		do
		{
			++index;
		} while (index < maximum_index_ && nb_refs(index) == 0);
		return index;
	}

protected:
	std::vector<AttributeGenT*> attributes_;
	std::vector<std::shared_ptr<AttributeGenT>> attributes_shared_ptr_;

	std::mutex mark_attributes_mutex_;
	std::vector<std::vector<AttributeGenT*>> mark_attributes_;
	std::vector<std::vector<uint32>> available_mark_attributes_;

	std::vector<uint32> available_indices_;

	uint32 nb_elements_;
	uint32 maximum_index_;

	friend AttributeGenT;

	void delete_attribute(AttributeGenT* attribute);

	virtual void init_ref_counter(uint32 index) = 0;
	virtual void reset_ref_counter(uint32 index) = 0;
	virtual uint32 nb_refs(uint32 index) const = 0;
	virtual void init_mark_attributes(uint32 index) = 0;
};

///////////////////////////////
// AttributeContainerT class //
///////////////////////////////

template <template <typename> class AttributeT>
class CGOGN_CORE_EXPORT AttributeContainerT : public AttributeContainerGen
{
public:
	template <typename T>
	using Attribute = AttributeT<T>;
	using AttributeGen = AttributeGenT;
	using MarkAttribute = Attribute<uint8>;

public:
	std::unique_ptr<Attribute<uint32>> ref_counter_;

	inline void init_ref_counter(uint32 index) override
	{
		static_cast<AttributeGenT*>(ref_counter_.get())
			->manage_index(index); // AttributeContainerT is friend of AttributeGenT
		(*ref_counter_)[index] = 1u;
	}

	inline void reset_ref_counter(uint32 index) override
	{
		(*ref_counter_)[index] = 0u;
	}

	inline uint32 nb_refs(uint32 index) const override
	{
		return (*ref_counter_)[index];
	}

	inline void init_mark_attributes(uint32 index) override
	{
		for (uint32 i = 0, nb = mark_attributes_.size(); i < nb; ++i)
		{
			for (AttributeGenT* mark_attribute : mark_attributes_[i])
			{
				MarkAttribute* m = static_cast<MarkAttribute*>(mark_attribute);
				(*m)[index] = 0u;
			}
		}
	}

public:
	AttributeContainerT() : AttributeContainerGen()
	{
		ref_counter_ = std::make_unique<Attribute<uint32>>(nullptr, "__refs");
	}

	~AttributeContainerT()
	{
	}

	template <typename T>
	std::shared_ptr<Attribute<T>> add_attribute(const std::string& name)
	{
		auto it = std::find_if(attributes_.begin(), attributes_.end(),
							   [&](AttributeGenT* att) { return att->name().compare(name) == 0; });
		if (it == attributes_.end())
		{
			std::shared_ptr<Attribute<T>> asp = std::make_shared<Attribute<T>>(this, name);
			Attribute<T>* ap = asp.get();
			static_cast<AttributeGenT*>(ap)->manage_index(
				maximum_index_); // AttributeContainerT is friend of AttributeGenT
			attributes_.push_back(ap);
			attributes_shared_ptr_.push_back(asp);
			return asp;
		}
		return std::shared_ptr<Attribute<T>>();
	}

	template <typename T>
	std::shared_ptr<Attribute<T>> get_attribute(const std::string& name) const
	{
		auto it = std::find_if(attributes_shared_ptr_.begin(), attributes_shared_ptr_.end(),
							   [&](const auto& att) { return att->name().compare(name) == 0; });
		if (it != attributes_shared_ptr_.end())
			return std::dynamic_pointer_cast<Attribute<T>>(*it);
		return std::shared_ptr<Attribute<T>>();
	}

	MarkAttribute* get_mark_attribute()
	{
		uint32 thread_index = current_thread_index();
		if (available_mark_attributes_[thread_index].size() > 0)
		{
			uint32 index = available_mark_attributes_[thread_index].back();
			available_mark_attributes_[thread_index].pop_back();
			return static_cast<MarkAttribute*>(mark_attributes_[thread_index][index]);
		}
		else
		{
			MarkAttribute* ap = new MarkAttribute(nullptr, "__mark");
			static_cast<AttributeGenT*>(ap)->manage_index(
				maximum_index_); // AttributeContainerT is friend of AttributeGenT
			mark_attributes_[thread_index].push_back(ap);
			return ap;
		}
	}

	void release_mark_attribute(MarkAttribute* attribute)
	{
		uint32 thread_index = current_thread_index();
		auto it = std::find(mark_attributes_[thread_index].begin(), mark_attributes_[thread_index].end(), attribute);
		cgogn_message_assert(it != mark_attributes_[thread_index].end(), "Mark Attribute not found on release");
		available_mark_attributes_[thread_index].push_back(std::distance(mark_attributes_[thread_index].begin(), it));
	}

	inline void ref_index(uint32 index)
	{
		cgogn_message_assert(nb_refs(index) > 0, "Trying to ref an unused index");
		(*ref_counter_)[index]++;
	}

	inline bool unref_index(uint32 index)
	{
		cgogn_message_assert(nb_refs(index) > 0, "Trying to unref an unused index");
		(*ref_counter_)[index]--;
		if ((*ref_counter_)[index] == 1u)
		{
			release_index(index);
			return true;
		}
		return false;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_CONTAINER_ATTRIBUTE_CONTAINER_H_

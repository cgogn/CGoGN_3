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

#ifndef CGOGN_CORE_CONTAINER_ATTRIBUTE_CONTAINER_H_
#define CGOGN_CORE_CONTAINER_ATTRIBUTE_CONTAINER_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/numerics.h>
#include <cgogn/core/utils/unique_ptr.h>

#include <vector>
#include <string>
#include <memory>

namespace cgogn
{

class AttributeContainerGen;
template <template <typename> class AttributeT>
class AttributeContainer;

////////////////////////
// AttributeGen class //
////////////////////////

class CGOGN_CORE_EXPORT AttributeGen
{
protected:

	AttributeContainerGen* container_;
	bool is_mark_;
	std::string name_;

public:

	AttributeGen(AttributeContainerGen* container, bool is_mark, const std::string& name);
	virtual ~AttributeGen();

	inline const std::string& name() const { return name_; }
	inline bool is_mark() const { return is_mark_; }

	uint32 maximum_index() const;

private:

	friend AttributeContainerGen;
	template <template <typename> class AttributeT> friend class AttributeContainer;

	virtual void manage_index(uint32 index) = 0;
};

/////////////////////////////////
// AttributeContainerGen class //
/////////////////////////////////

class CGOGN_CORE_EXPORT AttributeContainerGen
{
public:

	using AttributeGenPtr = std::shared_ptr<AttributeGen>;

protected:

	std::vector<AttributeGen*> attributes_;
	std::vector<AttributeGenPtr> attributes_shared_ptr_;

	std::vector<AttributeGen*> mark_attributes_;

	std::vector<uint32> available_indices_;

	uint32 nb_elements_;
	uint32 maximum_index_;

	friend AttributeGen;

	void delete_attribute(AttributeGen* attribute);

	virtual void init_ref_counter(uint32 index) = 0;
	virtual void reset_ref_counter(uint32 index) = 0;
	virtual uint32 nb_refs(uint32 index) const = 0;
	virtual void init_mark_attributes(uint32 index) = 0;

public:

	AttributeContainerGen();
	virtual ~AttributeContainerGen();

	inline uint32 nb_elements() const { return nb_elements_; }
	inline uint32 maximum_index() const { return maximum_index_; }

	uint32 new_index();
	void release_index(uint32 index);

	void remove_attribute(AttributeGenPtr attribute);

	using const_iterator = std::vector<AttributeGen*>::const_iterator;
	inline const_iterator begin() const { return attributes_.begin(); }
	inline const_iterator end() const { return attributes_.end(); }

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
		do { ++index; }
		while (
			index < maximum_index_ &&
			nb_refs(index) == 0
		);
		return index;
	}
};

//////////////////////////////
// AttributeContainer class //
//////////////////////////////

template <template <typename> class AttributeT>
class AttributeContainer : public AttributeContainerGen
{
public:

	template <typename T>
	using Attribute = AttributeT<T>;
	template <typename T>
	using AttributePtr = std::shared_ptr<AttributeT<T>>;

protected:

	std::unique_ptr<Attribute<uint32>> ref_counters_;

	void init_ref_counter(uint32 index) override
	{
		static_cast<AttributeGen*>(ref_counters_.get())->manage_index(index);
		(*ref_counters_)[index] = 1u;
	}

	void reset_ref_counter(uint32 index) override
	{
		(*ref_counters_)[index] = 0u;
	}

	void init_mark_attributes(uint32 index) override
	{
		for (auto mark_attribute : mark_attributes_)
		{
			Attribute<uint8>* m = static_cast<Attribute<uint8>*>(mark_attribute);
			(*m)[index] = 0u;
		}
	}

public:

	AttributeContainer() : AttributeContainerGen()
	{
		ref_counters_ = cgogn::make_unique<Attribute<uint32>>(this, true, "__refs");
	}

	~AttributeContainer()
	{}

	template <typename T>
	AttributePtr<T> add_attribute(const std::string& name)
	{
		AttributePtr<T> asp = get_attribute<T>(name);
		if (!asp)
			asp = std::make_shared<Attribute<T>>(this, false, name);
		Attribute<T>* ap = asp.get();
		static_cast<AttributeGen*>(ap)->manage_index(maximum_index_);
		attributes_.push_back(ap);
		attributes_shared_ptr_.push_back(asp);
		return asp;
	}

	template <typename T>
	AttributePtr<T> get_attribute(const std::string& name) const
	{
		auto it = std::find_if(
			attributes_shared_ptr_.begin(),
			attributes_shared_ptr_.end(),
			[&] (const AttributeGenPtr& att) { return att->name().compare(name) == 0; }
		);
		if (it != attributes_shared_ptr_.end())
			return std::dynamic_pointer_cast<Attribute<T>>(*it);
		return AttributePtr<T>();
	}

	std::unique_ptr<Attribute<uint8>> add_mark_attribute()
	{
		Attribute<uint8>* ap = new Attribute<uint8>(this, true, "__mark");
		static_cast<AttributeGen*>(ap)->manage_index(maximum_index_);
		mark_attributes_.push_back(ap);
		return std::unique_ptr<Attribute<uint8>>(ap);
	}

	uint32 nb_refs(uint32 index) const override
	{
		return (*ref_counters_)[index];
	}

	void ref_index(uint32 index)
	{
		cgogn_message_assert(nb_refs(index) > 0, "Trying to ref an unused index");
		(*ref_counters_)[index]++;
	}

	bool unref_index(uint32 index)
	{
		cgogn_message_assert(nb_refs(index) > 0, "Trying to unref an unused index");
		(*ref_counters_)[index]--;
		if ((*ref_counters_)[index] == 1u)
		{
			release_index(index);
			return true;
		}
		return false;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_CONTAINER_ATTRIBUTE_CONTAINER_H_
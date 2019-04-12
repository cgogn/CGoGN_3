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

	uint32 nb_elements_;
	uint32 maximum_index_;

	friend AttributeGen;

	void delete_attribute(AttributeGen* attribute);

public:

	using const_iterator = typename std::vector<AttributeGenPtr>::const_iterator;
	inline const_iterator begin() const { return attributes_shared_ptr_.begin(); }
	inline const_iterator end() const { return attributes_shared_ptr_.end(); }

	AttributeContainerGen();
	~AttributeContainerGen();

	inline uint32 nb_elements() const { return nb_elements_; }
	inline uint32 maximum_index() const { return maximum_index_; }

	uint32 get_index();

	void remove_attribute(AttributeGenPtr attribute);
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

	AttributeContainer() : AttributeContainerGen()
	{}

	~AttributeContainer()
	{}

	template <typename T>
	AttributePtr<T> add_attribute(const std::string& name)
	{
		AttributePtr<T> casp = std::make_shared<Attribute<T>>(this, false, name);
		Attribute<T>* cap = casp.get();
		static_cast<AttributeGen*>(cap)->manage_index(maximum_index_);
		attributes_.push_back(cap);
		attributes_shared_ptr_.push_back(casp);
		return casp;
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

	Attribute<uint8>* add_mark_attribute()
	{
		Attribute<uint8>* ca = new Attribute<uint8>(this, true, "mark");
		static_cast<AttributeGen*>(ca)->manage_index(maximum_index_);
		mark_attributes_.push_back(ca);
		return ca;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_CONTAINER_ATTRIBUTE_CONTAINER_H_

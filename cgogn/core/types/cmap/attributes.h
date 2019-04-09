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

#ifndef CGOGN_CORE_CMAP_ATTRIBUTES_H_
#define CGOGN_CORE_CMAP_ATTRIBUTES_H_

#include <cgogn/core/cgogn_core_export.h>

#include <cgogn/core/utils/numerics.h>

#include <vector>
#include <string>

namespace cgogn
{

class AttributeContainer;

class CGOGN_CORE_EXPORT AttributeGen
{
protected:

	std::string name_;
	AttributeContainer* container_;
	bool is_mark_;

	friend class AttributeContainer;
	virtual void add_line() = 0;
	virtual const void* data_ptr() const = 0;

public:

	AttributeGen(AttributeContainer* container, bool is_mark, const std::string& name);
	virtual ~AttributeGen();

	const std::string& name() const { return name_; }
};

template <typename T>
class CGOGN_CORE_EXPORT Attribute : public AttributeGen
{
	std::vector<T> data_;

	friend class AttributeContainer;
	void add_line() override { data_.push_back(T()); }
	void resize(uint32 size) { data_.resize(size); }

public:

	using const_iterator = typename std::vector<T>::const_iterator;
	inline const_iterator begin() const { return data_.begin(); }
	inline const_iterator end() const { return data_.end(); }

	using iterator = typename std::vector<T>::iterator;
	inline iterator begin() { return data_.begin(); }
	inline iterator end() { return data_.end(); }

	Attribute(AttributeContainer* container, bool is_mark, const std::string& name) : AttributeGen(container, is_mark, name)
	{}

	~Attribute() override
	{}

	uint32 size() const { return data_.size(); }
	const void* data_ptr() const override { return &data_[0]; }

	inline T& operator[](uint32 index) { return data_[index]; }
	inline const T& operator[](uint32 index) const { return data_[index]; }

	inline void swap(Attribute<T>* attribute)
	{
		if (attribute->container_ == this->container_)
			data_.swap(attribute->data_);
	}
};

class CGOGN_CORE_EXPORT AttributeContainer
{
	std::vector<AttributeGen*> attributes_;
	std::vector<Attribute<uint8>*> mark_attributes_;
	uint32 size_;

	friend class AttributeGen;

public:

	using const_iterator = std::vector<AttributeGen*>::const_iterator;
	inline const_iterator begin() const { return attributes_.begin(); }
	inline const_iterator end() const { return attributes_.end(); }

	AttributeContainer() : size_(0)
	{}

	~AttributeContainer()
	{
		for (AttributeGen* ag : attributes_)
			delete ag;
	}

	uint32 size() const { return size_; }

	template <typename T>
	Attribute<T>* add_attribute(const std::string& name)
	{
		Attribute<T>* a = new Attribute<T>(this, false, name);
		a->resize(size_);
		attributes_.push_back(a);
		return a;
	}

	template <typename T>
	Attribute<T>* get_attribute(const std::string& name) const
	{
		const_iterator it = std::find_if(attributes_.begin(), attributes_.end(), [&] (AttributeGen* att) { return att->name().compare(name) == 0; });
		if (it != attributes_.end())
		{
			Attribute<T>* res = dynamic_cast<Attribute<T>*>(*it);
			return res;
		}
		return nullptr;
	}

	Attribute<uint8>* add_mark_attribute()
	{
		Attribute<uint8>* a = new Attribute<uint8>(this, true, "mark");
		a->resize(size_);
		mark_attributes_.push_back(a);
		return a;
	}

	uint32 add_line()
	{
		for (AttributeGen* ag : attributes_)
			ag->add_line();
		for (Attribute<uint8>* a : mark_attributes_)
			a->add_line();
		return size_++;
	}
};

} // namespace cgogn

#endif // CGOGN_CORE_CMAP_ATTRIBUTES_H_

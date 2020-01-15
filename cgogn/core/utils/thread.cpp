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

#include <cgogn/core/utils/thread.h>
#include <cgogn/core/utils/thread_pool.h>

namespace cgogn
{

CGOGN_TLS uint32 thread_index_;
CGOGN_TLS Buffers<uint32>* uint32_buffers_thread_ = nullptr;

CGOGN_TLS uint32 uint32_value_ = 0;
CGOGN_TLS float64 float64_value_ = 0.0;

CGOGN_CORE_EXPORT uint32& uint32_value()
{
	return uint32_value_;
}
CGOGN_CORE_EXPORT float64& float64_value()
{
	return float64_value_;
}

CGOGN_CORE_EXPORT void thread_start(uint32 index)
{
	thread_index_ = index;
	if (uint32_buffers_thread_ == nullptr)
		uint32_buffers_thread_ = new Buffers<uint32>();
}

CGOGN_CORE_EXPORT void thread_stop()
{
	delete uint32_buffers_thread_;
	uint32_buffers_thread_ = nullptr;
}

CGOGN_CORE_EXPORT uint32 max_nb_threads()
{
	return uint32(thread_pool()->max_nb_workers()) + 1 // account for the main thread
		   + 1										   // account for 1 external thread
		;
}

CGOGN_CORE_EXPORT uint32 current_thread_index()
{
	return thread_index_;
}

CGOGN_CORE_EXPORT uint32 current_worker_index()
{
	return thread_index_ - 1; // account for the main thread
}

CGOGN_CORE_EXPORT Buffers<uint32>* uint32_buffers()
{
	return uint32_buffers_thread_;
}

} // namespace cgogn

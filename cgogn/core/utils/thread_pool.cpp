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

#include <cgogn/core/utils/thread_pool.h>
#include <cgogn/core/utils/thread.h>

namespace cgogn
{

ThreadPool::ThreadPool() :
    stop_(false)
{
	uint32 nb_ww = std::thread::hardware_concurrency() - 1;
	nb_working_workers_ = nb_ww;

	for (uint32 i = 0u; i < nb_ww; ++i)
	{
		workers_.emplace_back([this, i] () -> void
		{
			thread_start(i + 1);
			for(;;)
			{
				while (i >= nb_working_workers_)
				{
					std::unique_lock<std::mutex> lock(running_mutex_);
					condition_running_.wait(lock);
				}

				std::unique_lock<std::mutex> lock(queue_mutex_);
				condition_.wait(
					lock,
					[this] { return stop_ || !tasks_.empty(); }
				);

				if (stop_ && tasks_.empty())
				{
					thread_stop();
					return;
				}

				if (i < nb_working_workers_)
				{
					PackagedTask task = std::move(tasks_.front());
					tasks_.pop();
					lock.unlock();
#if defined(_MSC_VER) && _MSC_VER < 1900
					(*task)();
#else
					task();
#endif
				}
				else
				{
					lock.unlock();
					condition_.notify_one();
				}
			}
		});
	}

	std::cout << "ThreadPool launched with " << nb_working_workers_ << " workers" << std::endl;
}

ThreadPool::~ThreadPool()
{
	nb_working_workers_ = uint32(workers_.size());
	condition_running_.notify_all();

	{
		std::unique_lock<std::mutex> lock(queue_mutex_);
		stop_ = true;
	}
#if !(defined(CGOGN_WIN_VER) && (CGOGN_WIN_VER <= 61))
	condition_running_.notify_all();
	condition_.notify_all();
#endif
	for (std::thread& worker : workers_)
		worker.join();
}

void ThreadPool::set_nb_workers(uint32 nb)
{
	if (nb == 0xffffffff)
		nb_working_workers_ = uint32(workers_.size());
	else
		nb_working_workers_ = std::min(uint32(workers_.size()), nb);

	condition_running_.notify_all();

	std::cout << "ThreadPool now using " << nb_working_workers_ << " workers" << std::endl;
}

ThreadPool* thread_pool()
{
	// thread safe according to http://stackoverflow.com/questions/8102125/is-local-static-variable-initialization-thread-safe-in-c11
	static ThreadPool pool;
	return &pool;
}

} // namespace cgogn

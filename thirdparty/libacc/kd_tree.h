/*
 * Copyright (C) 2015-2018, Nils Moehrle
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef ACC_KDTREE_HEADER
#define ACC_KDTREE_HEADER

#include <algorithm>
#include <atomic>
#include <limits>
#include <queue>
#include <stack>
#include <thread>

// #include <math/vector.h>
// #include <array>
#include <Eigen/Dense>

#include "defines.h"

ACC_NAMESPACE_BEGIN

template <unsigned K, typename IdxType = unsigned>
class KDTree
{
	using Vector = Eigen::Matrix<double, K, 1>;

public:
	static constexpr IdxType NAI = std::numeric_limits<IdxType>::max();

	typedef std::shared_ptr<KDTree> Ptr;
	typedef std::shared_ptr<KDTree const> ConstPtr;

private:
	std::vector<Vector> const vertices;
	struct Node
	{
		typedef IdxType ID;
		IdxType vertex_id;
		decltype(K) d;
		union {
			IdxType first;
			Node::ID left;
		};
		union {
			IdxType last;
			Node::ID right;
		};
	};

	std::atomic<IdxType> num_nodes;
	std::vector<Node> nodes;
	typename Node::ID create_node(decltype(K) d, IdxType first, IdxType last)
	{
		typename Node::ID node_id = num_nodes++;
		Node& node = nodes[node_id];
		node.first = first;
		node.last = last;
		node.vertex_id = NAI;
		node.d = d;
		return node_id;
	}

	std::pair<typename Node::ID, typename Node::ID> ssplit(typename Node::ID node_id, std::vector<IdxType>* indices);

	void split(typename Node::ID node_id, std::vector<IdxType>* indices, std::atomic<int>* num_threads);

public:
	template <class C>
	static C convert(KDTree const& kd_tree);

	static Ptr create(std::vector<Vector> const& vertices, int max_threads = std::thread::hardware_concurrency())
	{
		return std::make_shared<KDTree>(vertices, max_threads);
	}

	KDTree(std::vector<Vector> const& vertices, int max_threads = std::thread::hardware_concurrency());

	bool find_nn(Vector point, std::pair<IdxType, double>* nn,
				 double max_dist = std::numeric_limits<double>::infinity()) const;

	bool find_nns(Vector point, std::size_t n, std::vector<std::pair<IdxType, double>>* nns_ptr,
				  double max_dist = std::numeric_limits<double>::infinity()) const;

	inline const Vector& vertex(IdxType idx) const
	{
		return vertices[idx];
	}
};

template <unsigned K, typename IdxType>
constexpr IdxType KDTree<K, IdxType>::NAI;

template <unsigned K, typename IdxType>
KDTree<K, IdxType>::KDTree(std::vector<Vector> const& vertices, int max_threads) : vertices(vertices), num_nodes(0)
{
	std::size_t num_vertices = vertices.size();
	nodes.resize(num_vertices);

	std::vector<IdxType> indices(num_vertices);
	std::iota(indices.begin(), indices.end(), 0);

	std::atomic<int> num_threads(max_threads);
	split(create_node(0, 0, num_vertices), &indices, &num_threads);
}

template <unsigned K, typename IdxType>
void KDTree<K, IdxType>::split(typename Node::ID node_id, std::vector<IdxType>* indices, std::atomic<int>* num_threads)
{
	typename Node::ID left, right;
	if ((*num_threads -= 1) >= 1)
	{
		std::tie(left, right) = ssplit(node_id, indices);
		if (left != NAI && right != NAI)
		{
			std::thread other(&KDTree::split, this, left, indices, num_threads);
			split(right, indices, num_threads);
			other.join();
		}
		else
		{
			if (left != NAI)
				split(left, indices, num_threads);
			if (right != NAI)
				split(right, indices, num_threads);
		}
	}
	else
	{
		std::deque<typename Node::ID> queue;
		queue.push_back(node_id);
		while (!queue.empty())
		{
			typename Node::ID node_id2 = queue.front();
			queue.pop_front();
			std::tie(left, right) = ssplit(node_id2, indices);
			if (left != NAI)
				queue.push_back(left);
			if (right != NAI)
				queue.push_back(right);
		}
	}
	*num_threads += 1;
}

template <unsigned K, typename IdxType>
std::pair<typename KDTree<K, IdxType>::Node::ID, typename KDTree<K, IdxType>::Node::ID> KDTree<K, IdxType>::ssplit(
	typename Node::ID node_id, std::vector<IdxType>* indices)
{
	Node& node = nodes[node_id];
	IdxType mid = (node.last + node.first) / 2;
	decltype(K) d = node.d;
	std::nth_element(indices->begin() + node.first, indices->begin() + mid, indices->begin() + node.last,
					 [this, d](IdxType a, IdxType b) -> bool { return vertices[a][d] < vertices[b][d]; });
	d = (d + 1) % K;
	node.vertex_id = indices->at(mid);
	if (mid - node.first > 0)
	{
		node.left = create_node(d, node.first, mid);
	}
	else
	{
		node.left = NAI;
	}
	if (node.last - (mid + 1) > 0)
	{
		node.right = create_node(d, mid + 1, node.last);
	}
	else
	{
		node.right = NAI;
	}
	return std::make_pair(node.left, node.right);
}

template <unsigned K, typename IdxType>
bool KDTree<K, IdxType>::find_nn(Vector point, std::pair<IdxType, double>* nn_ptr, double max_dist) const
{
	std::vector<std::pair<IdxType, double>> nns;
	if (!find_nns(point, 1, &nns, max_dist))
		return false;

	if (nn_ptr != nullptr)
		*nn_ptr = nns[0];
	return true;
}

template <typename IdxType>
static bool compare(std::pair<IdxType, double> l, std::pair<IdxType, double> r)
{
	return l.second < r.second;
}

template <unsigned K, typename IdxType>
bool KDTree<K, IdxType>::find_nns(Vector vertex, std::size_t n, std::vector<std::pair<IdxType, double>>* nns_ptr,
								  double max_dist) const
{
	std::vector<std::pair<IdxType, double>> nns;
	nns.reserve(n);

	// TODO use square distances

	typename Node::ID node_id = 0;
	std::stack<typename Node::ID> s;
	bool down = true;

	while (true)
	{
		Node const& node = nodes[node_id];

		double diff = vertex[node.d] - vertices[node.vertex_id][node.d];
		if (down)
		{
			double dist = (vertex - vertices[node.vertex_id]).norm();
			if (dist <= max_dist)
			{
				if (nns.size() < n)
				{
					nns.emplace_back(node.vertex_id, dist);
				}
				else
				{
					typename std::vector<std::pair<IdxType, double>>::iterator it;
					it = std::max_element(nns.begin(), nns.end(), compare<IdxType>);
					*it = std::make_pair(node.vertex_id, dist);
				}

				if (nns.size() == n)
				{
					typename std::vector<std::pair<IdxType, double>>::iterator it;
					it = std::max_element(nns.begin(), nns.end(), compare<IdxType>);
					max_dist = it->second;
				}
			}

			if (node.left != NAI || node.right != NAI)
			{
				/* Inner node - traverse further down. */
				down = true;

				if (node.left != NAI && node.right != NAI)
				{
					s.push(node_id);
				}

				double diff = vertex[node.d] - vertices[node.vertex_id][node.d];

				typename Node::ID next = (diff < 0.0f) ? node.left : node.right;
				typename Node::ID other = (diff < 0.0f) ? node.right : node.left;

				node_id = (next != NAI) ? next : other;
			}
			else
			{
				/* Leaf - traverse up and search for next node. */
				down = false;
				node_id = NAI;
			}
		}
		else
		{
			if (std::abs(diff) < max_dist)
			{
				down = true;
				node_id = (diff < 0.0f) ? node.right : node.left;
			}
			else
			{
				down = false;
				node_id = NAI;
			}
		}

		if (node_id == NAI)
		{
			if (s.empty())
				break;
			node_id = s.top();
			s.pop();
		}
	}

	std::sort(nns.begin(), nns.end(), compare<IdxType>);

	bool success = nns.size() == n;

	if (nns_ptr != nullptr)
		nns_ptr->swap(nns);

	return success;
}

ACC_NAMESPACE_END

#endif /* ACC_KDTREE_HEADER */

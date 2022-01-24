/*
 * Copyright (C) 2015-2018, Nils Moehrle
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef ACC_BVHTREE_HEADER
#define ACC_BVHTREE_HEADER

#include <algorithm>
#include <array>
#include <atomic>
#include <deque>
#include <limits>
#include <numeric>
#include <stack>
#include <thread>
#include <vector>

#include <cassert>

#include "primitives.h"

#ifndef BVHTREE_NUM_BINS
#define BVHTREE_NUM_BINS 64
#endif

ACC_NAMESPACE_BEGIN

template <typename IdxType, typename Vec3fType>
class BVHTree
{
public:
	typedef std::shared_ptr<BVHTree> Ptr;
	typedef std::shared_ptr<const BVHTree> ConstPtr;

	typedef acc::Ray<Vec3fType> Ray;
	struct Hit
	{
		/* Parameter of the ray (distance of hit location). */
		double t;
		/* Index of the struck triangle. */
		IdxType idx;
		/* Barycentric coordinates of hit location w.r.t. the triangle. */
		Vec3fType bcoords;
	};

private:
	static constexpr IdxType NAI = std::numeric_limits<IdxType>::max();

	typedef acc::AABB<Vec3fType> AABB;
	typedef acc::Tri<Vec3fType> Tri;

	struct Node
	{
		typedef IdxType ID;
		IdxType first;
		IdxType last;
		ID left;
		ID right;
		AABB aabb;
	};

	struct Bin
	{
		IdxType n;
		AABB aabb;
	};

	std::vector<IdxType> indices;
	std::vector<Tri> tris;

	std::atomic<IdxType> num_nodes;
	std::vector<Node> nodes;
	typename Node::ID create_node(IdxType first, IdxType last)
	{
		typename Node::ID node_id = num_nodes++;
		Node& node = nodes[node_id];
		node.first = first;
		node.last = last;
		node.left = NAI;
		node.right = NAI;
		node.aabb.min = Vec3fType(inf, inf, inf);
		node.aabb.max = Vec3fType(-inf, -inf, -inf);
		return node_id;
	}

	std::pair<typename Node::ID, typename Node::ID> sbsplit(typename Node::ID node_id, std::vector<AABB> const& aabbs);
	std::pair<typename Node::ID, typename Node::ID> bsplit(typename Node::ID node_id, std::vector<AABB> const& aabbs);
	std::pair<typename Node::ID, typename Node::ID> ssplit(typename Node::ID node_id, std::vector<AABB> const& aabbs);
	void split(typename Node::ID, std::vector<AABB> const& aabbs, std::atomic<int>* num_threads);

	bool intersect(Ray const& ray, typename Node::ID node_id, Hit* hit) const;
	std::pair<IdxType, Vec3fType> closest_point(Vec3fType vertex, typename Node::ID node_id) const;

public:
	static Ptr create(std::vector<IdxType> const& faces, std::vector<Vec3fType> const& vertices,
					  int max_threads = std::thread::hardware_concurrency())
	{
		return std::make_shared<BVHTree>(faces, vertices, max_threads);
	}

	template <class C>
	static C convert(BVHTree const& bvh_tree);

	/* Constructs the BVH tree using the Surface Area Heuristic as
	 * published in
	 * "On fast Construction of SAH-based Bounding Volume Hierarchies"
	 * by Ingo Wald (IEEE Symposium on Interactive Ray Tracing 2007)
	 *
	 * The mesh should be given as triangle index list and
	 * a vector containing the 3D positions. */
	BVHTree(std::vector<IdxType> const& faces, std::vector<Vec3fType> const& vertices,
			int max_threads = std::thread::hardware_concurrency());

	Tri const& get_triangle(IdxType idx) const
	{
		return tris[idx];
	}

	bool intersect(Ray ray, Hit* hit_ptr = nullptr) const;
	bool closest_point(Vec3fType vertex, std::pair<IdxType, Vec3fType>* cp_ptr, double max_dist = inf) const;
	Vec3fType closest_point(Vec3fType vertex);
};

template <typename IdxType, typename Vec3fType>
constexpr IdxType BVHTree<IdxType, Vec3fType>::NAI;

template <typename IdxType, typename Vec3fType>
void BVHTree<IdxType, Vec3fType>::split(typename Node::ID node, std::vector<AABB> const& aabbs,
										std::atomic<int>* num_threads)
{

	typename Node::ID left, right;
	if ((*num_threads -= 1) >= 1)
	{
		std::tie(left, right) = sbsplit(node, aabbs);
		if (left != NAI && right != NAI)
		{
			std::thread other(&BVHTree::split, this, left, std::cref(aabbs), num_threads);
			split(right, aabbs, num_threads);
			other.join();
		}
	}
	else
	{
		std::deque<typename Node::ID> queue;
		queue.push_back(node);
		while (!queue.empty())
		{
			typename Node::ID node2 = queue.back();
			queue.pop_back();
			std::tie(left, right) = sbsplit(node2, aabbs);
			if (left != NAI && right != NAI)
			{
				queue.push_back(left);
				queue.push_back(right);
			}
		}
	}
	*num_threads += 1;
}

template <typename IdxType, typename Vec3fType>
std::pair<typename BVHTree<IdxType, Vec3fType>::Node::ID, typename BVHTree<IdxType, Vec3fType>::Node::ID> BVHTree<
	IdxType, Vec3fType>::sbsplit(typename Node::ID node_id, std::vector<AABB> const& aabbs)
{
#if BVHTREE_NUM_BINS
	Node const& node = nodes[node_id];
	IdxType n = node.last - node.first;
	if (n > BVHTREE_NUM_BINS)
	{
		return bsplit(node_id, aabbs);
	}
	else
	{
		return ssplit(node_id, aabbs);
	}
#else
	return ssplit(node_id, aabbs);
#endif
}
#if BVHTREE_NUM_BINS

template <typename IdxType, typename Vec3fType>
std::pair<typename BVHTree<IdxType, Vec3fType>::Node::ID, typename BVHTree<IdxType, Vec3fType>::Node::ID> BVHTree<
	IdxType, Vec3fType>::bsplit(typename Node::ID node_id, std::vector<AABB> const& aabbs)
{
	Node& node = nodes[node_id];
	IdxType n = node.last - node.first;

	std::array<Bin, BVHTREE_NUM_BINS> bins;
	std::array<AABB, BVHTREE_NUM_BINS> right_aabbs;
	std::vector<IdxType> bin(n);

	double min_cost = inf;
	std::pair<IdxType, char> split;
	for (char d = 0; d < 3; ++d)
	{
		double min = node.aabb.min[d];
		double max = node.aabb.max[d];
		for (Bin& bin : bins)
		{
			bin = {0, {Vec3fType(inf, inf, inf), Vec3fType(-inf, -inf, -inf)}};
		}
		for (std::size_t i = node.first; i < node.last; ++i)
		{
			AABB const& aabb = aabbs[indices[i]];
			char idx = ((mid(aabb, d) - min) / (max - min)) * (BVHTREE_NUM_BINS - 1);
			bins[idx].aabb += aabb;
			bins[idx].n += 1;
			bin[i - node.first] = idx;
		}

		right_aabbs[BVHTREE_NUM_BINS - 1] = bins[BVHTREE_NUM_BINS - 1].aabb;
		for (std::size_t i = BVHTREE_NUM_BINS - 1; i > 0; --i)
		{
			right_aabbs[i - 1] = bins[i - 1].aabb + right_aabbs[i];
		}

		AABB left_aabb = bins[0].aabb;
		std::size_t nl = bins[0].n;
		for (std::size_t idx = 1; idx < BVHTREE_NUM_BINS; ++idx)
		{
			std::size_t nr = n - nl;
			double cost = 1.2f + (surface_area(left_aabb) / surface_area(node.aabb) * nl +
								  surface_area(right_aabbs[idx]) / surface_area(node.aabb) * nr);
			if (cost <= min_cost)
			{
				min_cost = cost;
				split = std::make_pair(d, idx);
			}

			nl += bins[idx].n;
			left_aabb += bins[idx].aabb;
		}
	}

	if (min_cost >= n)
		return std::make_pair(NAI, NAI);

	char d;
	IdxType sidx;
	std::tie(d, sidx) = split;

	double min = node.aabb.min[d];
	double max = node.aabb.max[d];
	for (Bin& bin : bins)
	{
		bin = {0, {Vec3fType(inf, inf, inf), Vec3fType(-inf, -inf, -inf)}};
	}

	for (std::size_t i = node.first; i < node.last; ++i)
	{
		AABB const& aabb = aabbs[indices[i]];
		char idx = ((mid(aabb, d) - min) / (max - min)) * (BVHTREE_NUM_BINS - 1);
		bins[idx].aabb += aabb;
		bins[idx].n += 1;
		bin[i - node.first] = idx;
	}

	IdxType l = node.first;
	IdxType r = node.last - 1;
	while (l < r)
	{
		if (bin[l - node.first] < sidx)
		{
			l += 1;
			continue;
		}
		if (bin[r - node.first] >= sidx)
		{
			r -= 1;
			continue;
		}
		std::swap(bin[l - node.first], bin[r - node.first]);
		std::swap(indices[l], indices[r]);
	}
	assert(l == r);
	std::size_t m = bin[(l & r) - node.first] >= sidx ? (l & r) : (l & r) + 1;

	node.left = create_node(node.first, m);
	node.right = create_node(m, node.last);
	for (std::size_t idx = 0; idx < BVHTREE_NUM_BINS; ++idx)
	{
		if (idx < sidx)
		{
			nodes[node.left].aabb += bins[idx].aabb;
		}
		else
		{
			nodes[node.right].aabb += bins[idx].aabb;
		}
	}

	return std::make_pair(node.left, node.right);
}
#endif

template <typename IdxType, typename Vec3fType>
std::pair<typename BVHTree<IdxType, Vec3fType>::Node::ID, typename BVHTree<IdxType, Vec3fType>::Node::ID> BVHTree<
	IdxType, Vec3fType>::ssplit(typename Node::ID node_id, std::vector<AABB> const& aabbs)
{
	Node& node = nodes[node_id];
	IdxType n = node.last - node.first;

	double min_cost = inf;
	std::pair<char, IdxType> split;
	std::vector<AABB> right_aabbs(n);
	for (char d = 0; d < 3; ++d)
	{
		std::sort(&indices[node.first], &indices[node.last], [&aabbs, d](IdxType first, IdxType second) -> bool {
			return mid(aabbs[first], d) < mid(aabbs[second], d) ||
				   (mid(aabbs[first], d) == mid(aabbs[second], d) && first < second);
		});

		right_aabbs[n - 1] = aabbs[indices[node.last - 1]];
		for (IdxType i = node.last - 1; i > node.first; --i)
		{
			right_aabbs[i - 1 - node.first] = aabbs[indices[i - 1]] + right_aabbs[i - node.first];
		}
		node.aabb = right_aabbs[0];

		AABB left_aabb = aabbs[indices[node.first]];
		for (IdxType i = node.first + 1; i < node.last; ++i)
		{
			IdxType nl = i - node.first;
			IdxType nr = n - nl;
			double cost = 1.2f + (surface_area(left_aabb) / surface_area(node.aabb) * nl +
								  surface_area(right_aabbs[nl]) / surface_area(node.aabb) * nr);
			if (cost <= min_cost)
			{
				min_cost = cost;
				split = std::make_pair(d, i);
			}

			left_aabb += aabbs[indices[i]];
		}
	}

	if (min_cost >= n)
		return std::make_pair(NAI, NAI);

	char d;
	IdxType i;
	std::tie(d, i) = split;
	std::sort(&indices[node.first], &indices[node.last], [&aabbs, d](std::size_t first, std::size_t second) -> bool {
		return mid(aabbs[first], d) < mid(aabbs[second], d) ||
			   (mid(aabbs[first], d) == mid(aabbs[second], d) && first < second);
	});

	node.left = create_node(node.first, i);
	node.right = create_node(i, node.last);
	return std::make_pair(node.left, node.right);
}

template <typename IdxType, typename Vec3fType>
BVHTree<IdxType, Vec3fType>::BVHTree(std::vector<IdxType> const& faces, std::vector<Vec3fType> const& vertices,
									 int max_threads)
	: num_nodes(0)
{

	std::size_t num_faces = faces.size() / 3;
	std::vector<AABB> aabbs(num_faces);
	std::vector<Tri> ttris(num_faces);

	/* Initialize vector with upper bound of nodes. */
	nodes.resize(2 * num_faces - 1);

	/* Initialize root node. */
	Node& root = nodes[create_node(0, num_faces)];
	for (std::size_t i = 0; i < aabbs.size(); ++i)
	{
		ttris[i].a = vertices[faces[i * 3 + 0]];
		ttris[i].b = vertices[faces[i * 3 + 1]];
		ttris[i].c = vertices[faces[i * 3 + 2]];

		calculate_aabb(ttris[i], &aabbs[i]);
		root.aabb += aabbs[i];
	}

	indices.resize(num_faces);
	std::iota(indices.begin(), indices.end(), 0);

	std::atomic<int> num_threads(max_threads);
	split(0, aabbs, &num_threads);

	nodes.resize(num_nodes);

	/* Order breath first */
	typename Node::ID node_id = 0;
	typename Node::ID idx = 0;

	IdxType last = 0;

	std::vector<IdxType> nindices;
	nindices.reserve(indices.size());
	std::vector<Node> nnodes(nodes.size());

	nnodes[idx++] = nodes[node_id];
	while (node_id < nnodes.size())
	{
		Node& node = nnodes[node_id++];
		if (node.left != NAI && node.right != NAI)
		{
			nnodes[idx] = nodes[node.left];
			node.left = idx++;
			nnodes[idx] = nodes[node.right];
			node.right = idx++;
		}
		else
		{
			nindices.insert(nindices.begin() + last, indices.begin() + node.first, indices.begin() + node.last);
			node.first = last;
			last = static_cast<IdxType>(nindices.size());
			node.last = last;
		}
	}

	std::swap(nnodes, nodes);
	std::swap(nindices, indices);

	tris.resize(ttris.size());
	for (std::size_t i = 0; i < indices.size(); ++i)
		tris[i] = ttris[indices[i]];
}

template <typename IdxType, typename Vec3fType>
bool BVHTree<IdxType, Vec3fType>::intersect(Ray const& ray, typename Node::ID node_id, Hit* hit) const
{
	Node const& node = nodes[node_id];
	bool ret = false;
	for (std::size_t i = node.first; i < node.last; ++i)
	{
		double t;
		Vec3fType bcoords;
		if (acc::intersect(ray, tris[i], &t, &bcoords))
		{
			if (t > hit->t)
				continue;
			hit->idx = indices[i];
			hit->t = t;
			hit->bcoords = bcoords;
			ret = true;
		}
	}
	return ret;
}

template <typename IdxType, typename Vec3fType>
bool BVHTree<IdxType, Vec3fType>::intersect(Ray ray, Hit* hit_ptr) const
{
	Hit hit;
	hit.t = inf;

	typename Node::ID node_id = 0;
	std::stack<typename Node::ID> s;
	while (true)
	{
		Node const& node = nodes[node_id];
		if (node.left != NAI && node.right != NAI)
		{
			double tmin_left, tmin_right;
			bool left = acc::intersect(ray, nodes[node.left].aabb, &tmin_left);
			bool right = acc::intersect(ray, nodes[node.right].aabb, &tmin_right);
			if (left && right)
			{
				if (tmin_left < tmin_right)
				{
					s.push(node.right);
					node_id = node.left;
				}
				else
				{
					s.push(node.left);
					node_id = node.right;
				}
			}
			else
			{
				if (right)
					node_id = node.right;
				if (left)
					node_id = node.left;
			}

			if (!left && !right)
			{
				if (s.empty())
					break;
				node_id = s.top();
				s.pop();
			}
		}
		else
		{
			if (intersect(ray, node_id, &hit))
				ray.tmax = hit.t;

			if (s.empty())
				break;
			node_id = s.top();
			s.pop();
		}
	}

	if (hit.t < inf)
	{
		if (hit_ptr != nullptr)
			*hit_ptr = hit;
		return true;
	}
	else
		return false;
}

template <typename IdxType, typename Vec3fType>
std::pair<IdxType, Vec3fType> BVHTree<IdxType, Vec3fType>::closest_point(Vec3fType vertex,
																		 typename Node::ID node_id) const
{
	Node const& node = nodes[node_id];

	IdxType idx = NAI;
	Vec3fType cp;
	double dist = inf;

	for (std::size_t i = node.first; i < node.last; ++i)
	{
		Vec3fType cp_tri = acc::closest_point(vertex, tris[i]);
		double dist_tri = (cp_tri - vertex).squaredNorm();
		if (dist_tri < dist)
		{
			cp = cp_tri;
			dist = dist_tri;
			idx = i;
		}
	}

	return std::make_pair(idx, cp);
}

template <typename IdxType, typename Vec3fType>
bool BVHTree<IdxType, Vec3fType>::closest_point(Vec3fType vertex, std::pair<IdxType, Vec3fType>* cp_ptr,
												double max_dist) const
{
	double dist = max_dist * max_dist;
	IdxType idx = NAI;
	Vec3fType cp;

	typename Node::ID node_id = 0;
	std::stack<typename Node::ID> s;
	while (true)
	{
		Node const& node = nodes[node_id];
		if (node.left != NAI && node.right != NAI)
		{
			Vec3fType cp_left = acc::closest_point(vertex, nodes[node.left].aabb);
			Vec3fType cp_right = acc::closest_point(vertex, nodes[node.right].aabb);
			double dmin_left = (cp_left - vertex).squaredNorm();
			double dmin_right = (cp_right - vertex).squaredNorm();
			bool left = dmin_left < dist;
			bool right = dmin_right < dist;
			if (left && right)
			{
				if (dmin_left < dmin_right)
				{
					s.push(node.right);
					node_id = node.left;
				}
				else
				{
					s.push(node.left);
					node_id = node.right;
				}
			}
			else
			{
				if (right)
					node_id = node.right;
				if (left)
					node_id = node.left;
			}

			if (!left && !right)
			{
				if (s.empty())
					break;
				node_id = s.top();
				s.pop();
			}
		}
		else
		{
			IdxType idx_leaf;
			Vec3fType cp_leaf;
			std::tie(idx_leaf, cp_leaf) = closest_point(vertex, node_id);
			double dist_leaf = (cp_leaf - vertex).squaredNorm();
			if (dist_leaf < dist)
			{
				dist = dist_leaf;
				idx = idx_leaf;
				cp = cp_leaf;
			}

			if (s.empty())
				break;
			node_id = s.top();
			s.pop();
		}
	}

	if (idx == NAI)
		return false;

	if (cp_ptr != nullptr)
	{
		cp_ptr->first = indices[idx];
		cp_ptr->second = cp;
	}

	return true;
}

template <typename IdxType, typename Vec3fType>
Vec3fType BVHTree<IdxType, Vec3fType>::closest_point(Vec3fType vertex)
{
	std::pair<IdxType, Vec3fType> cp;
	closest_point(vertex, &cp);
	return cp.second;
}

ACC_NAMESPACE_END

#endif /* ACC_BVHTREE_HEADER */

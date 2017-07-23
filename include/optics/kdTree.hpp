// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)


#pragma once

#include <array>
#include <algorithm>
#include <fplus/fplus.hpp>



namespace kdt{


constexpr bool is_powerof2(std::size_t v) {
	return v && ((v & (v - 1)) == 0);
}


template<typename CoordsType, std::size_t dimension, std::size_t n_points, std::size_t max_points_per_node, std::size_t split_dim = 0 >
class KDTree {
	//static_assert(is_powerof2(n_points), "The total number of points stored in a KDTree 'n_points' has to be a power of two.");
	//static_assert(is_powerof2(max_pts_per_node), "");

public:
	typedef std::array<CoordsType, dimension> point_type;

	KDTree() {}

	KDTree(const std::array<point_type, n_points>& points_)
	{
		//std::array<std::size_t, n_points> indices = fplus::numbers<std::size_t, std::array<std::size_t, n_points>>(static_cast<std::size_t>(0), n_points);
		//fplus::nth_element 
		/*auto split_index_it = std::begin(indices);
		std::advance(split_index_it, static_cast<typename split_index_it::difference_type>(n));
		auto split_index_it = std::nth_element(std::begin(indices), split_index_it, std::end(indices), comp<split_dim>);
		auto split_index = *split_index_it;
		auto split_point = points[split_index];
		split = split_point[split_dim];
		*/
		const auto comp = [](const point_type& p1, const point_type& p2) -> bool {
			return p1[split_dim] < p2[split_dim];
		};

		auto split_point = fplus::nth_element_by(comp, n_points / 2, points_);
		split = split_point[split_dim];

		std::array<point_type, n_points / 2> left_points; //Could hold references instead of real points?
		std::array<point_type, n_points - n_points / 2> right_points;
		std::size_t l_idx(0);
		std::size_t r_idx(0);
		for (const auto& p : points_) {
			if (p[split_dim] < split) {
				left_points[l_idx++] = p;
			}
			else {
				right_points[r_idx++] = p;
			}
		}
		if (n_points/2 <= max_points_per_node) {
			left = std::make_unique<KDTreeLeaf<CoordsType, dimension, n_points / 2, n_points / 2, (split_dim + 1) % dimension>>(left_points);
			right = std::make_unique<KDTreeLeaf<CoordsType, dimension, n_points - n_points / 2, n_points - n_points / 2, (split_dim + 1) % dimension>>( right_points );
		}
		else {
			left = std::make_unique<KDTree<CoordsType, dimension, n_points / 2, max_points_per_node, (split_dim + 1) % dimension>>(left_points);
			right = std::make_unique<KDTree<CoordsType, dimension, n_points - n_points / 2, max_points_per_node, (split_dim + 1) % dimension>> (right_points);
		
		}
	}

	virtual std::vector<std::array<point_type, max_points_per_node>> radius_search (const point_type& p, double radius) {
		if (p[split_dim] + radius < split) {
			return left->radius_search(p, radius);
		}
		else if (p[split_dim] - radius >= split) {
			return right->radius_search(p, radius);
		}
		else {
			return fplus::append(left->radius_search(p, radius), right->radius_search(p, radius));
		}
	}


private:
	double split;
	std::unique_ptr<KDTreeLeaf<CoordsType, dimension, n_points/2, max_points_per_node, (split_dim + 1) % dimension>> left;
	std::unique_ptr<KDTreeLeaf<CoordsType, dimension, n_points - (n_points / 2), max_points_per_node, (split_dim + 1) % dimension>> right;
};


template<typename CoordsType, std::size_t dimension, std::size_t n_points, std::size_t max_points_per_node, std::size_t split_dim >
class KDTreeLeaf : public KDTree<CoordsType, dimension, n_points, max_points_per_node, split_dim> {
public:
	KDTreeLeaf(const std::array<point_type, max_points_per_node>& points_) : points(points_) {}

	std::vector<std::array<point_type, max_points_per_node>> radius_search (const point_type& p, double radius) override {
		return std::vector<std::array<point_type, max_points_per_node>>({ points });
	}

private:
	std::array<point_type, max_points_per_node> points;
};

}//namespace kdt

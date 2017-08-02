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
class KDTree;

template<typename CoordsType, std::size_t dimension>
static double square_distance( const std::array<CoordsType, dimension>& p1, const std::array<CoordsType, dimension>& p2 ) {
	double result = 0.0;
	for ( std::size_t i = 0; i < dimension; i++ ) {
		double d = (p1[i] - p2[i]);
		result +=  d*d;
	}
	return result;
}


template<typename CoordsType, std::size_t dimension, std::size_t n_points >
class KDTreeLeaf {

	/*
	template<CoordsType, dimension, 2 * n_point, std::size_t max_points_per_node, std::size_t split_dim>
	friend class KDTree;
	template<CoordsType, dimension, 2 * n_points + 1, std::size_t max_points_per_node, std::size_t split_dim>
	friend class KDTree;
	*/
public:
	typedef std::array<CoordsType, dimension> Point;
	KDTreeLeaf() {}
	KDTreeLeaf( const std::array<Point, n_points>& points_ ) : points( points_ ) {}

	void radius_search_( const Point& p, double radius, std::vector<Point>& neighbors ){//TODO: make private & friend KDTree
		for ( const auto& x : points ) {
			if ( square_distance( p, x ) <= radius*radius ) {
				neighbors.push_back( x );
			}
		}
	}

private:
	std::array<Point, n_points> points;
};


/*
template<typename CoordsType, std::size_t dimension, std::size_t n_points, std::size_t max_points_per_node, std::size_t split_dim>
std::enable_if_t<
	(n_points <= max_points_per_node),
	std::unique_ptr<KDTree<CoordsType, dimension, n_points, max_points_per_node, (split_dim + 1) % dimension>>
>
make_child( const std::array<std::array<CoordsType, dimension>, n_points>& points_ ) {
	return std::make_unique<KDTreeLeaf<CoordsType, dimension, n_points, max_points_per_node, (split_dim + 1) % dimension>>( points_ );
}

template<typename CoordsType, std::size_t dimension, std::size_t n_points, std::size_t max_points_per_node, std::size_t split_dim>
std::enable_if_t<
	(n_points > max_points_per_node),
	std::unique_ptr<KDTree<CoordsType, dimension, n_points, max_points_per_node, (split_dim + 1) % dimension>>
>
make_child( const std::array<std::array<CoordsType, dimension>, n_points>& points_ ) {
	return std::make_unique<KDTree<CoordsType, dimension, n_points, max_points_per_node, (split_dim + 1) % dimension>>( points_ );
}
*/

template<typename CoordsType, std::size_t dimension, std::size_t n_points, std::size_t max_points_per_node, std::size_t split_dim>
class KDTree {
	//static_assert(is_powerof2(n_points), "The total number of points stored in a KDTree 'n_points' has to be a power of two.");
	//static_assert(is_powerof2(max_pts_per_node), "");

public:
	typedef std::array<CoordsType, dimension> Point;

	KDTree() {}

	KDTree(const std::array<Point, n_points>& points_)
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
		const auto comp = []( const Point& p1, const Point& p2 ) {
			return p1[split_dim] < p2[split_dim];
		};

		auto split_point = fplus::nth_element_by(comp, n_points / 2, points_);
		split = split_point[split_dim];

		std::array<Point, n_points / 2> left_points; //Could hold references instead of real points?
		std::array<Point, n_points - n_points / 2> right_points;
		std::size_t l_idx(0);
		std::size_t r_idx(0);
		std::vector<Point> on_split_points;
		on_split_points.reserve( 2 );
		for (const auto& p : points_) {
			if (l_idx < left_points.size() && p[split_dim] < split) {
				left_points[l_idx++] = p;
			}
			else if ( p[split_dim] > split ){
				right_points[r_idx++] = p;
			}
			else{
				on_split_points.push_back( p );
			}
		}
		assert( left_points.size() + right_points.size() == l_idx + r_idx + on_split_points.size() );

		for ( const auto& p : on_split_points ) {
			if ( l_idx < left_points.size() ) {
				left_points[l_idx++] = p;
			}
			else {
				right_points[r_idx++] = p;
			}
		}
		
		left = decltype(left)( left_points );
		right = decltype(right)( right_points );
	
	}


	std::vector<Point> radius_search (const Point& p, double radius) {
		std::vector<Point> neighbors;
		neighbors.reserve( 2* max_points_per_node ); //TODO: intelligent guess of number of expected neighbors?
		radius_search_( p, radius, neighbors );
		neighbors.shrink_to_fit();
		return neighbors;
	}


//private:
	double split;

//	constexpr std::size_t n_left_points() { return n_points / 2; }
//	constexpr std::size_t n_right_points() { return n_points - (n_points / 2) };

	using LeftChild =
		typename std::conditional<(n_points / 2 <= max_points_per_node),
		KDTreeLeaf<CoordsType, dimension, n_points / 2>,
		KDTree<CoordsType, dimension, n_points / 2, max_points_per_node, (split_dim + 1) % dimension>
		>::type;
	using RightChild =
		typename std::conditional< (n_points - (n_points / 2) <= max_points_per_node),
		KDTreeLeaf<CoordsType, dimension, n_points - (n_points / 2)>,
		KDTree<CoordsType, dimension, n_points - (n_points / 2), max_points_per_node, (split_dim + 1) % dimension>
		>::type;

	//<CoordsType, dimension, n_points / 2, max_points_per_node, (split_dim + 1) % dimension>
	//<CoordsType, dimension, n_points / 2, max_points_per_node, (split_dim + 1) % dimension>
	LeftChild left;
	RightChild right;

	void radius_search_( const Point& p, double radius, std::vector<Point>& neighbors ) { //TODO: make private
		if ( p[split_dim] + radius < split ) {
			return left.radius_search_( p, radius, neighbors );
		}
		else if ( p[split_dim] - radius >= split ) {
			return right.radius_search_( p, radius, neighbors );
		}
		else {
			left.radius_search_( p, radius, neighbors );
			right.radius_search_( p, radius, neighbors );
			return;
		}
	}

};



}//namespace kdt

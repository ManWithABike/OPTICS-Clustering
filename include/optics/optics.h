// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)


#pragma once

#define _HAS_AUTO_PTR_ETC 1

#include "bgr_image.hpp"

#include <geometry/geometry.hpp>
#include <fplus/fplus.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include<vector>
#include <exception>



namespace optics {


struct reachability_dist {
	reachability_dist( std::size_t point_index, double reach_dist ) : point_index( point_index ), reach_dist( reach_dist ) {}

    std::string to_string() const{
        return std::to_string( point_index) + ":" + std::to_string( reach_dist );
    }
	std::size_t point_index;
	double reach_dist;
};

inline bool operator<( const reachability_dist& lhs, const reachability_dist& rhs ) {
	return (lhs.reach_dist == rhs.reach_dist) ? (lhs.point_index < rhs.point_index) : (lhs.reach_dist < rhs.reach_dist);
}
inline bool operator==( const reachability_dist& lhs, const reachability_dist& rhs ) {
	return (lhs.reach_dist == rhs.reach_dist) && (lhs.point_index == rhs.point_index) ;
}


namespace internal{


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

template<typename T, std::size_t N>
using Pt = typename bg::model::point<T, N, bg::cs::cartesian>;

template<typename T, std::size_t N>
using TreeValue = typename std::pair<Pt<T, N>, std::size_t>;

template<typename T, std::size_t N>
using Box = typename bg::model::box<Pt<T, N>>;

template<typename T, std::size_t N>
using RTree = typename bgi::rtree<TreeValue<T, N>, bgi::quadratic<25>>; //TODO: Number of elems per node configurable?



//recursive set boost coordinates template
template <typename T, size_t N, size_t I>
struct set_boost_point_coords
{
	static inline int set( Pt<T, N>& boost_pt, const geom::Vec<T, N>& coords )
	{
		bg::set<I>( boost_pt, coords[I] );
		return set_boost_point_coords<T, N, I - 1>::set( boost_pt, coords );
	}
};

template <typename T, size_t N>
struct set_boost_point_coords<T, N, 0>
{
	static inline int set( Pt<T, N>& boost_pt, const geom::Vec<T, N>& coords )
	{
		bg::set<0>( boost_pt, coords[0] );
		return 0;
	}
};

template<typename T, std::size_t N>
Pt<T, N> geom_to_boost_point( const geom::Vec<T, N>& point ) {
	Pt<T, N> boost_point;
	set_boost_point_coords<T, N, N - 1 >::set( boost_point, point );
	return boost_point;
}



template<typename T, std::size_t N>
RTree<T, N> initialize_rtree( const std::vector<geom::Vec<T, N>>& points ) {
	//Insert all points with index into cloud
	size_t idx_id = 0;
	auto cloud = fplus::transform( [&idx_id]( const geom::Vec<T, N>& point ) -> TreeValue<T, N> {
		return{ geom_to_boost_point( point ),  idx_id++ };
	}, points );

	//Create an rtree from the cloud using packaging
	RTree<T, N> rtree( cloud );
	return rtree;
}

template<typename T, std::size_t dimension>
double dist( const Pt<T,dimension>& boost_pt, const geom::Vec<T, dimension>& geom_pt ) {
	const auto dist = bg::distance( boost_pt, geom_to_boost_point( geom_pt ) );
	return dist;
}

template<typename T, std::size_t N>
std::vector<std::size_t> find_neighbor_indices( const geom::Vec<T, N>& point, const double epsilon, const RTree<T, N>& rtree ) {
	//produce search box
	geom::Vec<double, N> eps_vec( epsilon );
	geom::Vec<double, N> corner1 = point.as_doubles() - eps_vec;
	geom::Vec<double, N> corner2 = point.as_doubles() + eps_vec;
	Box<double, N> query_box( geom_to_boost_point( corner1 ), geom_to_boost_point( corner2 ) ); //TODO: possible speed-up if query-box is of type T (e.g. int) when possible?

	//search neighbors in box (manhattan dist 2*epsilon)
	std::vector<TreeValue<T, N>> neighbors;
	rtree.query( bgi::intersects( query_box ), std::back_inserter( neighbors ) );

	//keep those with euclidean dist < epsilon
	auto neighbor_indices = fplus::transform_and_keep_justs( [&point, &epsilon]( const TreeValue<T, N>& tree_val ) ->  fplus::maybe<std::size_t> {
		if ( dist<T,N>( tree_val.first, point ) > epsilon ) {
			return {};
		}
		return tree_val.second;
	}, neighbors );

	return neighbor_indices;
}


template<typename T, std::size_t N>
fplus::maybe<double> compute_core_dist( const geom::Vec<T, N>& point, const std::vector<geom::Vec<T, N>>& points, const std::vector<std::size_t>& neighbor_indices, const std::size_t min_pts ) {

	if ( neighbor_indices.size() < min_pts ) { return{}; }

	//sort neighbors by distance
	/*auto neighbors = fplus::sort_on( [&points, &point]( const std::size_t& idx ) -> double {
		return geom::square_dist( point, points[idx] );
	}, neighbor_indices );

	double core_dist = geom::dist( point, points[neighbors[min_pts - 1]] );
	*/
	auto core_elem_idx = fplus::nth_element_on( [&points, &point]( const std::size_t& idx ) -> double {
		return geom::square_dist( point, points[idx] );
	}, min_pts-1, neighbor_indices );
	double core_dist = geom::dist( points[core_elem_idx], point );
	return core_dist;
}


inline void erase_idx_from_set( const reachability_dist& d, std::set<reachability_dist>& seeds ) {
	auto x = seeds.erase( d );
	assert( x == 1 );
}

template<typename T>
T pop_from_set( std::set<T>& set ) {
	T element = *set.begin();
	set.erase( set.begin() );
	return element;
}


template<typename T, std::size_t N>
void update( const geom::Vec<T, N>& point, const std::vector<geom::Vec<T, N>>& points, const std::vector<std::size_t>& neighbor_indices, const double core_dist,
				const std::vector<bool>& processed, std::vector<double>& reachability, std::set<reachability_dist>& seeds
			) {
	for ( const auto& o : neighbor_indices ) {
		if ( processed[o] ) { continue; }
		double new_reachability_dist = fplus::max( core_dist, geom::dist( point, points[o] ) );
		if ( reachability[o] < 0.0 ) {
			reachability[o] = new_reachability_dist;
			seeds.insert( reachability_dist( o, new_reachability_dist ) );
		}

		else if ( new_reachability_dist < reachability[o] ) {
			//erase from seeds
			erase_idx_from_set( reachability_dist( o, reachability[o]), seeds );
			//update reachability
			reachability[o] = new_reachability_dist;
			//reinsert seed with new reachability
			seeds.insert( reachability_dist( o, new_reachability_dist ) );
		}
	}
}

} //namespace internal



template<typename T, std::size_t dimension>
std::vector<reachability_dist> compute_reachability_dists( const std::vector<geom::Vec<T, dimension>>& points, const std::size_t min_pts, double epsilon ) {
	static_assert(std::is_convertible<double,T>::value, "optics::compute_reachability_dists: Point type 'T' must be convertible to double!" );
	static_assert( dimension >= 1, "optics::compute_reachability_dists: dimension must be >=1");
	if ( points.empty() ) { return{}; }

	//algorithm tracker
	std::vector<bool> processed( points.size(), false );
	std::vector<std::size_t> ordered_list;
	ordered_list.reserve( points.size() );
	std::vector<double> reachability( points.size(), -1.0f );
	//std::vector<double> core_dist( points.size(), -1.0f );

	//the rtree for fast nearest neighbour search
	auto rtree = internal::initialize_rtree( points );

	for ( std::size_t point_idx = 0; point_idx < points.size(); point_idx++ ) {
		if ( processed[point_idx] == true ) continue;
		processed[point_idx] = true;
		ordered_list.push_back( point_idx );
		std::set<reachability_dist> seeds;

		auto neighbor_indices = internal::find_neighbor_indices( points[point_idx], epsilon, rtree );

		fplus::maybe<double> core_dist_m = internal::compute_core_dist( points[point_idx], points, neighbor_indices, min_pts );
		if ( !core_dist_m.is_just() ) { continue; }
		double core_dist = core_dist_m.unsafe_get_just();

		internal::update( points[point_idx], points, neighbor_indices, core_dist, processed, reachability, seeds );

		while ( !seeds.empty() ) {
			reachability_dist s = internal::pop_from_set( seeds );
			assert( processed[s.point_index] == false );
			processed[s.point_index] = true;
			ordered_list.push_back( s.point_index );

			auto s_neighbor_indices = internal::find_neighbor_indices( points[s.point_index], epsilon, rtree );

			auto s_core_dist_m = internal::compute_core_dist( points[s.point_index], points, s_neighbor_indices, min_pts );
			if ( !s_core_dist_m.is_just() ) { continue; }
			double s_core_dist = s_core_dist_m.unsafe_get_just();

			internal::update( points[s.point_index], points, s_neighbor_indices, s_core_dist, processed, reachability, seeds );
		}

	}

	//sanity checks
	assert( ordered_list.size() == points.size() );
	assert( fplus::all_unique( ordered_list ) );

	//merge reachabilities into ordered list
	auto result = fplus::transform( [&reachability]( std::size_t point_idx ) -> reachability_dist {
		return reachability_dist( point_idx, reachability[point_idx] );
	}, ordered_list );
	return result;
}


template<typename T, std::size_t dimension>
std::vector<reachability_dist> compute_reachability_dists( const std::vector<std::array<T, dimension>>& points, const std::size_t min_pts, double epsilon ) {
	static_assert(std::is_convertible<double, T>::value, "optics::compute_reachability_dists: Point type 'T' must be convertible to double!");
	static_assert(dimension >= 1, "optics::compute_reachability_dists: dimension must be >=1");
	if ( points.empty() ) { return{}; }

	std::vector<geom::Vec<T, dimension>> geom_points;
	geom_points.reserve( points.size() );
	for ( const auto& p : points ) {
		geom_points.push_back( geom::Vec<T, dimension>( p ) );
	}

	return compute_reachability_dists( geom_points, min_pts, epsilon );
}


inline std::vector<std::vector<std::size_t>> get_cluster_indices( const std::vector<reachability_dist>& reach_dists, double reachability_threshold ) {
	assert( reach_dists.front().reach_dist < 0.0 );
	std::vector<std::vector<std::size_t>> result;
	for ( const auto& r : reach_dists ) {
		if ( r.reach_dist < 0.0 || r.reach_dist >= reachability_threshold ) {
			result.push_back( { r.point_index } );
		}
		else {
			result.back().push_back( r.point_index );
		}
	}
	return result;
}

template<typename T, std::size_t dimension>
std::vector<std::vector<std::array<T, dimension>>> get_cluster_points( const std::vector<reachability_dist>& reach_dists, double reachability_threshold, const std::vector<std::array<T, dimension>>& points ) {
	const auto sorted_reachdist_indices = fplus::unique( fplus::sort_on( []( const reachability_dist& r )->std::size_t { return r.point_index; }, reach_dists ) );
	assert( sorted_reachdist_indices.size() == points.size() );
	assert( sorted_reachdist_indices.back().point_index == points.size() - 1 );

	auto clusters = get_cluster_indices( reach_dists, reachability_threshold );
	std::vector<std::vector<std::array<T, dimension>>> result;
	result.reserve( clusters.size() );
	for ( const auto& cluster_indices : clusters ) {
		result.push_back( fplus::elems_at_idxs( cluster_indices, points ) );
	}
	return result;
}

template<typename T, std::size_t dimension>
std::vector<std::vector<geom::Vec<T, dimension>>> get_cluster_points( const std::vector<reachability_dist>& reach_dists, double reachability_threshold, const std::vector<geom::Vec<T, dimension>>& points ) {
	const auto sorted_reachdists = fplus::unique( fplus::sort( reach_dists ) );
	assert( sorted_reachdists.size() == points.size() );
	assert( sorted_reachdists.back() == points.size() - 1 );

	auto clusters = get_cluster_indices( reach_dists, reachability_threshold );
	std::vector<geom::Vec<T, dimension>> result;
	result.reserve( clusters.size() );
	for ( const auto& cluster_indices : clusters ) {
		result.push_back( fplus::elems_at_idxs( cluster_indices, points ) );
	}
	return result;
}

/**Exports a list of reachability distances into a csv file.
* If switch replace_nodists is set (as by default), points that haven't been assigned a reachability distance will be set to the maximum reachability distance + 1.
* Otherwise, points that haven't been assigned a reachability distance at all will appear with reachability distance -1.0.
*/
inline void export_reachability_dists( const std::vector<reachability_dist>& reach_dists, const std::string& csv_file_path, bool replace_nodists = true ) {
	//open filestream
	std::ofstream stream;
	stream.open( csv_file_path, std::ios::binary );
	if ( !stream.is_open() ) {
		std::cout << "export_reachability_dists(): Stream to file \"" + csv_file_path + "\" could not be opened!";
		return;
	}

	double no_dist = -1;
	if ( replace_nodists ) {
		no_dist = fplus::maximum_on( []( const reachability_dist& r ) -> double {
			return r.reach_dist;
		}, reach_dists ).reach_dist;
		no_dist += 1;
	}

	//write csv header
	stream << "PointIndex;ReachabilityDistance" << std::endl;
	//write body
	for ( const auto& x : reach_dists ) {
		stream << x.point_index << ";" << (x.reach_dist < 0 ? no_dist : x.reach_dist) << std::endl;
	}
}



inline void draw_reachability_plot( const std::vector<reachability_dist>& reach_dists, const std::string& img_file_name ) {
	if ( reach_dists.size() < 2 ) return;
	bgr_image image( size_2d( std::max( reach_dists.size(), std::size_t( 100 ) ), 256 ), bgr_col(255,255,255) );

	auto reach_dists_values = fplus::transform( []( const reachability_dist& r )-> double {
		return r.reach_dist;
	}, reach_dists );

	//Normalize the data
	double max_val = fplus::maximum( reach_dists_values );
	reach_dists_values.push_back( max_val + fplus::max( 30, max_val / 3 ) );//The future no_dist for points which weren't assigned any reachability dist. Has to be at least 30, and scale with max_val. Will be normalized to 256-64
	reach_dists_values.push_back( 1.0 );//In order to see where 1.0 was mapped after the normalization
	reach_dists_values = fplus::normalize_min_max( -1, 256 - 64, reach_dists_values );
	//Extract normalized 1.0:
	double one = reach_dists_values.back();
	reach_dists_values.pop_back();
	//Extract no_dist:
	int no_dist = fplus::min( 255, fplus::round(reach_dists_values.back()));
	reach_dists_values.pop_back();

	for ( int i = 0; i < static_cast<int>(reach_dists.size())-1; i++ ) {
        int x1 = fplus::round( (image.size().width_-1) * i/ static_cast<double>((reach_dists_values.size()-1)) );
        int y1 = image.size().height_ -1 - (reach_dists_values[i] < 0 ? no_dist : fplus::round( reach_dists_values[i]));
        int x2 = fplus::round( (image.size().width_-1) * (i+1)/static_cast<double>((reach_dists_values.size()-1)) );
        int y2 = image.size().height_ -1 - (reach_dists_values[i+1] < 0 ? no_dist : fplus::round( reach_dists_values[i+1]));
		plot_line_segment( image, img_pos(x1, y1), img_pos(x2, y2), bgr_col(0,0,0) );
		bgr_col col = reach_dists_values[i] < 0 ? bgr_col(0, 0, 255) : bgr_col(0, 255, 0);
		draw_pixel( image, img_pos( x1, y1 ), col );
		col = reach_dists_values[i+1] < 0 ? bgr_col( 0, 0, 255 ) : bgr_col( 0, 255, 0 );
		draw_pixel( image, img_pos( x2, y2 ), col );
	}

	//Draw Scale
	int x2 = fplus::round( image.size().width_ / static_cast<double>((reach_dists.size() - 1)) );
	plot_line_segment( image, img_pos( 0, image.size().height_ - 1 ), img_pos( x2, image.size().height_ - 1 ), bgr_col(0,255,0) ); //One point
	plot_line_segment( image, img_pos( 0, image.size().height_ - 1 ), img_pos( 0, image.size().height_ - 1 - fplus::round<double,std::size_t>( one ) ), bgr_col( 255, 0, 0 ) ); //ReachDist 1
	int no_dist_marker = fplus::min( static_cast<int>(image.size().height_ -1), fplus::round( fplus::maximum( reach_dists_values ) + 10.0 ));
	plot_line_segment( image, img_pos( 0, image.size().height_ - 1 - no_dist_marker ), img_pos( image.size().width_/3, image.size().height_ - 1 - no_dist_marker ), bgr_col( 0, 0, 255 ) ); //ReachDist 1

    image.save(img_file_name);
	return;
}

template<typename T>
bgr_image draw_2d_clusters( const std::vector<std::vector<geom::Vec<T,2>>>& clusters ) {
	auto box = geom2d::bounding_box( fplus::concat( clusters ) );
	bgr_image cluster_image( size_2d(box.get_size().first+1, box.get_size().second+1), bgr_col(255,255,255) );
	std::array<bgr_col, 6> colours = { bgr_col( 255,0,0 ), bgr_col( 0,255,0 ), bgr_col(0,0,255), bgr_col(255,255,0), bgr_col(255,0,255), bgr_col(0,255,255) };
	int col_idx = 0;
	for ( const auto& cluster : clusters ) {
		bgr_col col = colours[col_idx];
		++col_idx %= colours.size();
		auto cluster_box = geom2d::bounding_box( cluster );
		for( const auto& edge: fplus::overlapping_pairs_cyclic(cluster_box.points())){
            plot_line_segment( cluster_image, img_pos(edge.first.x(), edge.first.y()), img_pos(edge.second.x(), edge.second.y()), col );
		}
		for ( const auto & pt : cluster ) {
			cluster_image.pix( img_pos( pt.x()-box.bl().x(), pt.y()-box.bl().y() ) ) = col;
		}
	}
	return cluster_image;
}

template<typename T>
bgr_image draw_2d_clusters( const std::vector<std::vector<std::array<T, 2>>>& clusters ) {
	auto geom_clusters = fplus::transform_inner( []( const std::array<T, 2>& pt )-> geom::Vec<T, 2> {
		return geom::Vec<T,2>( pt );
	}, clusters );

	return draw_2d_clusters( geom_clusters );
}

} //namespace optics

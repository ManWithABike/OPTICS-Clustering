// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)

#include "CImg/CImg.h"
#include "geometry/geometry.hpp"

#include <fplus/fplus.hpp>

#pragma warning( push )
#pragma warning( disable: 4503 )
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#pragma warning( pop )

#include<vector>
#include <exception>



namespace optics {


struct reachability_dist {
	reachability_dist( std::size_t point_index, double reach_dist ) : point_index( point_index ), reach_dist( reach_dist ) {}

    std::string to_string() const{
        return std::to_string( point_index) + std::to_string( reach_dist );
    }
	std::size_t point_index;
	double reach_dist;
};

inline bool operator<( const reachability_dist& lhs, const reachability_dist& rhs ) {
	assert( lhs.point_index != rhs.point_index );
	return (lhs.reach_dist == rhs.reach_dist) ? (lhs.point_index < rhs.point_index) : (lhs.reach_dist < rhs.reach_dist);
}



namespace internal{


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

template<typename T, std::size_t N>
using Pt = typename bg::model::point<T, N, bg::cs::cartesian>;

template<typename T, std::size_t N>
using TreeValues = typename std::pair<Pt<T, N>, std::size_t>;

template<typename T, std::size_t N>
using Box = typename bg::model::box<Pt<T, N>>;

template<typename T, std::size_t N>
using RTree = typename bgi::rtree<TreeValues<T, N>, bgi::quadratic<25>>; //TODO: Number of elems per node configurable?



//recursive set boost coordinates template
template <typename T, size_t N, size_t I>
struct set_boost_point_coords
{
	static inline int set( Pt<T, N>& boost_pt, const geom::Vec<T, N>& coords )
	{
		bg::set<I>( boost_pt, coords[I] );
		return 0;
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
	set_boost_point_coords<T, N, N - 1>::set( boost_point, point );
	return boost_point;
}


template<typename T, std::size_t N>
RTree<T, N> initialize_rtree( const std::vector<geom::Vec<T, N>>& points ) {
	//Insert all points with index into cloud
	size_t idx_id = 0;
	auto cloud = fplus::transform( [&idx_id]( const geom::Vec<T, N>& point ) -> TreeValues<T, N> {
		return{ geom_to_boost_point( point ),  idx_id++ };
	}, points );

	//Create an rtree from the cloud using packaging
	RTree<T, N> rtree( cloud );
	return rtree;
}


template<typename T, std::size_t N>
std::vector<std::size_t> find_neighbor_indices( const geom::Vec<T, N>& point, std::vector<geom::Vec<T, N>> points, T epsilon, const RTree<T, N>& rtree ) {
	auto point_boost = geom_to_boost_point( point );
	//produce search box
	geom::Vec<T, N> eps_vec( epsilon );
	geom::Vec<T, N> corner1 = point - eps_vec;
	geom::Vec<T, N> corner2 = point + eps_vec;
	Box<T, N> query_box( geom_to_boost_point( corner1 ), geom_to_boost_point( corner2 ) );

	//search neighbors in box (manhattan dist 2*epsilon)
	std::vector<TreeValues<T, N>> neighbors;
	rtree.query( bgi::intersects( query_box ), std::back_inserter( neighbors ) );

	//extract indices from tree values
	std::vector<std::size_t> neighbor_indices = fplus::transform( []( const TreeValues<T, N>& v ) -> std::size_t {
		return v.second;
	}, neighbors );

	//keep those with euclidean dist < epsilon
	std::vector<std::size_t> drop_indices;
	for ( std::size_t i = 0; i < neighbor_indices.size(); i++ ) {
		if ( geom::dist( point, points[neighbor_indices[i]] ) > epsilon ) {
			drop_indices.push_back( i );
		}
	}
	neighbor_indices = fplus::drop_idxs( drop_indices, neighbor_indices );
	return neighbor_indices;
}


template<typename T, std::size_t N>
fplus::maybe<double> compute_core_dist( const geom::Vec<T, N>& point, const std::vector<geom::Vec<T, N>>& points, const std::vector<std::size_t>& neighbor_indices, const double epsilon, const std::size_t min_pts ) {

	if ( neighbor_indices.size() < min_pts ) { return{}; }

	//sort neighbors by distance
	auto neighbors = fplus::sort_on( [&points, &point]( const std::size_t& idx ) -> double {
		return geom::dist( point, points[idx] );
	}, neighbor_indices );

	double core_dist = geom::dist( point, points[neighbors[min_pts - 1]] );
	return core_dist;
}


inline void erase_idx_from_set( const std::size_t idx, std::set<reachability_dist>& seeds ) {
	for ( auto it = seeds.begin(); it != seeds.end(); it++ ) {
		if ( it->point_index == idx ) {
			seeds.erase( it );
			return;
		}
	}
	assert( false );
}

template<typename T>
T pop_from_set( std::set<T>& set ) {
	T element = *set.begin();
	set.erase( set.begin() );
	return element;
}


template<typename T, std::size_t N>
void update( const geom::Vec<T, N>& point, const std::vector<geom::Vec<T, N>>& points, const std::vector<std::size_t>& neighbor_indices, const double core_dist,
				const std::vector<bool>& processed, std::vector<double>& reachability, std::set<reachability_dist>& seeds,
				const double epsilon, const std::size_t min_pts ) {
	for ( const auto& o : neighbor_indices ) {
		if ( processed[o] ) { continue; }
		double new_reachability_dist = fplus::max( core_dist, geom::dist( point, points[o] ) );
		if ( reachability[o] < 0.0 ) {
			reachability[o] = new_reachability_dist;
			seeds.insert( reachability_dist( o, new_reachability_dist ) );
		}

		else if ( new_reachability_dist < reachability[o] ) {
			reachability[o] = new_reachability_dist;
			//loesche idx aus seeds
			erase_idx_from_set( o, seeds );
			//fuege idx mit new_reach_dist wieder in seeds ein
			seeds.insert( reachability_dist( o, new_reachability_dist ) );
		}
	}
}

} //namespace internal



template<typename T, std::size_t N>
std::vector<reachability_dist> compute_reachability_dists( const std::vector<geom::Vec<T, N>>& points, const std::size_t min_pts, double epsilon ) {
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

		auto neighbor_indices = internal::find_neighbor_indices( points[point_idx], points, epsilon, rtree );

		fplus::maybe<double> core_dist_m = internal::compute_core_dist( points[point_idx], points, neighbor_indices, epsilon, min_pts );
		if ( !core_dist_m.is_just() ) { continue; }
		double core_dist = core_dist_m.unsafe_get_just();

		internal::update( points[point_idx], points, neighbor_indices, core_dist, processed, reachability, seeds, epsilon, min_pts );

		while ( !seeds.empty() ) {
			reachability_dist s = internal::pop_from_set( seeds );
			assert( processed[s.point_index] == false );
			processed[s.point_index] = true;
			ordered_list.push_back( s.point_index );

			auto s_neighbor_indices = internal::find_neighbor_indices( points[s.point_index], points, epsilon, rtree );

			auto s_core_dist_m = internal::compute_core_dist( points[s.point_index], points, s_neighbor_indices, epsilon, min_pts );
			if ( !s_core_dist_m.is_just() ) { continue; }
			double s_core_dist = s_core_dist_m.unsafe_get_just();

			internal::update( points[s.point_index], points, s_neighbor_indices, s_core_dist, processed, reachability, seeds, epsilon, min_pts );
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
	if ( points.empty() ) { return{}; }
	assert( dimension != 0 );

	std::vector<geom::Vec<T, dimension>> geom_points(points.size());
	for ( const auto& p : points ) {
		if ( p.size() != dimension ) {
			assert(false); //TODO:Exception
			//throw(std::exception( std::string("compute_reachability_dists(): All Points must have the same dimension (i.e. number of coordinates) in order to be clustered.") ));
		}
		geom_points.push_back( geom::Vec<T, dimension>( p ) );
	}

	return compute_reachability_dists( geom_points, min_pts, epsilon );
}


inline std::vector<std::vector<std::size_t>> make_clusters( const std::vector<reachability_dist>& reach_dists, double reachability_threshold ) {
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


/**Exports a list of reachability distances into a csv file.
* If switch replace_nodists is set (as by default), points that haven't been assigned a reachability distance will be set to the maximum reachability distance + 1.
* Otherwise, points that haven't been assigned a reachability distance at all will appear with reachability distance -1.0.
*/
inline void export_reachability_dists( const std::vector<reachability_dist>& reach_dists, const std::string& csv_file_path, double replace_nodists = true ) {
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
	namespace cil = cimg_library;

	double no_dist = fplus::maximum_on( []( const reachability_dist& r ) -> double {
		return r.reach_dist;
	}, reach_dists ).reach_dist + 1;

	cil::CImg<double> graph_img( reach_dists.size(), 1, 1, 1 );
	for ( std::size_t i = 0; i < reach_dists.size(); i++ ) {
		graph_img( i, 0, 0, 0 ) = reach_dists[i].reach_dist < 0 ? no_dist : reach_dists[i].reach_dist;
	}

	cil::CImg<double> plot_img( std::max( reach_dists.size(), std::size_t(100)), 128, 1, 1, 0 );
	float col = 255.0;
	float* col_ptr = &col;

	plot_img.draw_graph( graph_img, col_ptr );

	plot_img.save( img_file_name.c_str() );
	return;
}

} //namespace optics

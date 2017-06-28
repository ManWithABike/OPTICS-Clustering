// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)


#pragma once

#ifndef _HAS_AUTO_PTR_ETC
#define _HAS_AUTO_PTR_ETC 1
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#undef _HAS_AUTO_PTR_ETC
#else
static_assert(_HAS_AUTO_PTR_ETC, "_HAS_AUTO_PTR_ETC has to be 1 for boost includes in MSVC_17, but has externally already been set to 0");
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#endif

#include "bgr_image.hpp"
#include "Tree.hpp"

#include <geometry/geometry.hpp>
#include <fplus/fplus.hpp>

#include <vector>
#include <exception>




namespace optics {


struct reachability_dist {
	reachability_dist( std::size_t point_index_, double reach_dist_ ) : point_index( point_index_ ), reach_dist( reach_dist_ ) {}

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
double epsilon_estimation( const std::vector<geom::Vec<T, dimension>>& points, const std::size_t min_pts ){
	static_assert(std::is_convertible<double, T>::value, "optics::epsilon_estimation: Point type 'T' must be convertible to double!");
	static_assert(dimension >= 1, "optics::epsilon_estimation: dimension must be >=1");
	if ( points.empty() ) { return 0; }

	double d = static_cast<double> (dimension);
	auto  space = geom::bounding_box( points );
	double space_volume = geom::product( geom::abs(space.first - space.second) );
	double nominator = space_volume * static_cast<double>(min_pts) * std::tgamma( d/2.0 + 1.0 );
	double denominator = static_cast<double>(points.size()) * std::sqrt( std::pow( geom::pi, d ) );
	double r = std::pow( nominator / denominator, 1.0 / d );
	return r;
}


template<typename T, std::size_t dimension>
double epsilon_estimation( const std::vector<std::array<T, dimension>>& points, const std::size_t min_pts ) {
	static_assert(std::is_convertible<double, T>::value, "optics::epsilon_estimation: Point type 'T' must be convertible to double!");
	static_assert(dimension >= 1, "optics::epsilon_estimation: dimension must be >=1");
	if ( points.empty() ) { return 0; }

	std::vector<geom::Vec<T, dimension>> geom_points;
	geom_points.reserve( points.size() );
	for ( const auto& p : points ) {
		geom_points.push_back( geom::Vec<T, dimension>( p ) );
	}

	return epsilon_estimation( geom_points, min_pts );
}


template<typename T, std::size_t dimension>
std::vector<reachability_dist> compute_reachability_dists( const std::vector<geom::Vec<T, dimension>>& points, const std::size_t min_pts, double epsilon = 0.0 ) {
	static_assert(std::is_convertible<T,double>::value, "optics::compute_reachability_dists: Point type 'T' must be convertible to double!" );
	static_assert( dimension >= 1, "optics::compute_reachability_dists: dimension must be >=1");
	if ( points.empty() ) { return{}; }

	if ( epsilon <= 0.0 ) {
		epsilon = epsilon_estimation( points, min_pts );
	}

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
std::vector<reachability_dist> compute_reachability_dists( const std::vector<std::array<T, dimension>>& points, const std::size_t min_pts, double epsilon = 0.0 ) {
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


bgr_image draw_reachability_plot( const std::vector<reachability_dist>& reach_dists, std::size_t min_width = 100 ) {
	if ( reach_dists.size() < 2 ) return bgr_image(size_2d(0,0));
	bgr_image image( size_2d( std::max( reach_dists.size(), min_width ), 256 ), bgr_col(255,255,255) );

	auto reach_dists_values = fplus::transform( []( const reachability_dist& r )-> double {
		return r.reach_dist;
	}, reach_dists );

	//Normalize the data
	double max_val = fplus::maximum( reach_dists_values );
	reach_dists_values.push_back( max_val + fplus::max( 30, max_val / 3 ) );//The future no_dist for points which weren't assigned any reachability dist. Has to be at least 30, and scale with max_val. Will be normalized to 256-64
	reach_dists_values.push_back( 10.0 );//In order to see where 10.0 was mapped after the normalization
	reach_dists_values.push_back( -1.0 );//In order to have at least one no dist 
	reach_dists_values = fplus::normalize_min_max( -1.0, 256.0 - 64.0, reach_dists_values );
	reach_dists_values.pop_back();
	//Extract normalized 10.0:
	double ten = reach_dists_values.back();
	reach_dists_values.pop_back();
	//Extract no_dist:
	int no_dist = fplus::min( 255, fplus::round(reach_dists_values.back()));
	reach_dists_values.pop_back();

	//Drawing the graph
	for ( int i = 0; i < static_cast<int>(reach_dists.size())-1; i++ ) {
        int x1 = fplus::round( (image.size().width_-1) * i/ static_cast<double>((reach_dists_values.size()-1)) );
        int y1 = image.size().height_ -1 - (reach_dists_values[i] < 0 ? no_dist : fplus::round( reach_dists_values[i]));
        int x2 = fplus::round( (image.size().width_-1) * (i+1)/static_cast<double>((reach_dists_values.size()-1)) );
        int y2 = image.size().height_ -1 - (reach_dists_values[i+1] < 0 ? no_dist : fplus::round( reach_dists_values[i+1]));
		plot_line_segment( image, img_pos(x1, y1), img_pos(x2, y2), bgr_col(30,30,30) );
		bgr_col col = reach_dists_values[i] < 0 ? bgr_col(0, 0, 255) : bgr_col(0, 255, 0);
		plot_pixel( image, img_pos( x1, y1 ), col );
		col = reach_dists_values[i+1] < 0 ? bgr_col( 0, 0, 255 ) : bgr_col( 0, 255, 0 );
		plot_pixel( image, img_pos( x2, y2 ), col );
	}

	//Fill area under the graph
	bgr_col fill_col( 177, 177, 177 );
	for ( std::size_t x = 0; x < image.size().width_; x++ ) {
		int y = image.size().height_-1; 
		while ( y >= 0 && image.pix( img_pos( x, y ) ) == bgr_col( 255, 255, 255 ) ) {
			image.pix( img_pos( x, y ) ) = fill_col;
			y--;
		}
	}

	//Draw Scale
	int x2 = fplus::round( image.size().width_ / static_cast<double>((reach_dists.size() - 1)) );
	plot_line_segment( image, img_pos( 0, image.size().height_ - 1 ), img_pos( x2, image.size().height_ - 1 ), bgr_col(0,255,0) ); //One point
	plot_line_segment( image, img_pos( 0, image.size().height_ - 1 ), img_pos( 0, image.size().height_ - 1 - fplus::round<double,std::size_t>( ten ) ), bgr_col( 255, 0, 0 ) ); //ReachDist 1
	int no_dist_marker = fplus::min( static_cast<int>(image.size().height_ -1), fplus::round( fplus::maximum( reach_dists_values ) + 10.0 ));
	plot_line_segment( image, img_pos( 0, image.size().height_ - 1 - no_dist_marker ), img_pos( image.size().width_/3, image.size().height_ - 1 - no_dist_marker ), bgr_col( 0, 0, 255 ) ); //ReachDist 1

	return image;
}


inline std::vector<std::vector<std::size_t>> get_cluster_indices( const std::vector<reachability_dist>& reach_dists, double reachability_threshold, const std::string& reachability_plot_image_path = "") {
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
    if( reachability_plot_image_path!= ""){
        auto img = draw_reachability_plot( reach_dists);
		img.save( reachability_plot_image_path );
	}
	return result;
}

template<typename T, std::size_t dimension>
std::vector<std::vector<std::array<T, dimension>>> get_cluster_points(
            const std::vector<reachability_dist>& reach_dists,
            double reachability_threshold,
            const std::vector<std::array<T, dimension>>& points,
            const std::string& reachability_plot_image_path = "" )
{
	const auto sorted_reachdist_indices =
        fplus::unique(
            fplus::sort_on(
                []( const reachability_dist& r )->std::size_t { return r.point_index; },
                reach_dists
            )
        );
	assert( sorted_reachdist_indices.size() == points.size() );
	assert( sorted_reachdist_indices.back().point_index == points.size() - 1 );

	auto clusters = get_cluster_indices( reach_dists, reachability_threshold );
	std::vector<std::vector<std::array<T, dimension>>> result;
	result.reserve( clusters.size() );
	for ( const auto& cluster_indices : clusters ) {
		result.push_back( fplus::elems_at_idxs( cluster_indices, points ) );
	}

	if( reachability_plot_image_path!= ""){
        auto img = draw_reachability_plot( reach_dists);
		img.save( reachability_plot_image_path );
	}
	return result;
}


template<typename T, std::size_t dimension>
std::vector<std::vector<geom::Vec<T, dimension>>> get_cluster_points(
        const std::vector<reachability_dist>& reach_dists,
        double reachability_threshold,
        const std::vector<geom::Vec<T, dimension>>& points ,
        const std::string& reachability_plot_image_path = "" )
{
	const auto sorted_reachdists = fplus::unique( fplus::sort( reach_dists ) );
	assert( sorted_reachdists.size() == points.size() );
	assert( sorted_reachdists.back() == points.size() - 1 );

	auto clusters = get_cluster_indices( reach_dists, reachability_threshold );
	std::vector<geom::Vec<T, dimension>> result;
	result.reserve( clusters.size() );
	for ( const auto& cluster_indices : clusters ) {
		result.push_back( fplus::elems_at_idxs( cluster_indices, points ) );
	}
	if( reachability_plot_image_path!= ""){
        auto img = draw_reachability_plot( reach_dists );
		img.save( reachability_plot_image_path );
	}
	return result;
}


struct SDA{
    SDA( std::size_t begin_idx, std::size_t end_idx, double mib) : begin_idx(begin_idx), end_idx(end_idx), mib(mib){}
    std::size_t begin_idx;
    std::size_t end_idx;
    double mib;
};


typedef std::pair<std::size_t, std::size_t> chi_cluster_indices;
typedef optics::Tree<chi_cluster_indices> cluster_tree;


std::vector<chi_cluster_indices> get_chi_clusters_flat( const std::vector<reachability_dist>& reach_dists, const double chi, std::size_t min_pts ){
    std::vector<std::pair<std::size_t, std::size_t>> clusters;
    std::vector<SDA> SDAs;
    double mib(0);
    double max_reach(0.0);
	for( const auto& r :reach_dists ){ if( r.reach_dist > max_reach ) max_reach = r.reach_dist; }
    const auto get_reach_dist = [&reach_dists, &max_reach](const std::size_t idx) -> double{
        assert( idx <= reach_dists.size() );
        if( idx == reach_dists.size() ) return max_reach;
		if( idx == 0 ) return max_reach;
        return reach_dists[idx].reach_dist;
    };
    const auto is_steep_down_pt = [&reach_dists, &chi](std::size_t idx ){
        if( idx == 0 ) return true;
        if( idx+1 >= reach_dists.size() ) return false;
        return reach_dists[idx+1].reach_dist <= reach_dists[idx].reach_dist * (1-chi);
    };
    const auto is_steep_up_pt = [&reach_dists, &chi](std::size_t idx ){
        if( idx+1 >= reach_dists.size() ) return true;
        return reach_dists[idx+1].reach_dist * (1-chi) >= reach_dists[idx].reach_dist;
    };
    const auto filter_sdas = [&chi, &reach_dists, &SDAs, &mib, &get_reach_dist](){
        SDAs = fplus::keep_if( [&mib, &reach_dists, &chi, &get_reach_dist](const SDA& sda)->bool{
                    return mib <= get_reach_dist(sda.begin_idx) * (1-chi);
            }, SDAs);
		for ( auto& sda : SDAs ) {
			sda.mib = std::max( sda.mib, mib );
		}
    };
    const auto get_sda_end = [&chi, &reach_dists, &min_pts, &is_steep_down_pt]( const std::size_t start_idx ) -> std::size_t{
        assert( is_steep_down_pt(start_idx) );
        std::size_t last_sd_idx = start_idx;
        std::size_t idx = start_idx +1;
        while( idx < reach_dists.size() ){
            if( idx - last_sd_idx >= min_pts){ return last_sd_idx; }
            if( reach_dists[idx].reach_dist > reach_dists[idx-1].reach_dist ){ return last_sd_idx; }
            if( is_steep_down_pt(idx) ){ last_sd_idx = idx; }
            idx++;
        }
		return std::max( reach_dists.size() - 2, last_sd_idx);
	};
	const auto get_sua_end = [&chi, &reach_dists, &min_pts, &is_steep_up_pt]( const std::size_t start_idx ) -> std::size_t {
		assert( is_steep_up_pt( start_idx ) );
		std::size_t last_su_idx = start_idx;
		std::size_t idx = start_idx + 1;
		while ( idx < reach_dists.size() ) {
			if ( idx - last_su_idx >= min_pts ) { return last_su_idx; }
			if ( reach_dists[idx].reach_dist < reach_dists[idx - 1].reach_dist ) { return last_su_idx; }
			if ( is_steep_up_pt( idx ) ) { last_su_idx = idx; }
			idx++;
		}
		return std::max( reach_dists.size() - 2, last_su_idx );
	};
	const auto cluster_borders = [&reach_dists, &chi]( const SDA& sda, std::size_t sua_begin_idx, std::size_t sua_end_idx ) -> std::pair<std::size_t, std::size_t> {
		double start_reach = reach_dists[sda.begin_idx].reach_dist;
		double end_reach = reach_dists[std::min( sua_end_idx + 1, reach_dists.size() - 1 )].reach_dist;
		if ( geom::in_range( start_reach, end_reach, start_reach*chi ) ) {
			return { sda.begin_idx, sua_end_idx };
		}
		if ( start_reach > end_reach ) {
			std::size_t start_idx = sda.begin_idx + 1;
			while ( start_idx <= sda.end_idx && reach_dists[start_idx].reach_dist > end_reach ) {
				start_idx++;
			}
			return { start_idx - 1, sua_end_idx };
		}
		if ( start_reach < end_reach ) {
			std::size_t end_idx = sua_end_idx;
			while ( end_idx >= sua_begin_idx && reach_dists[end_idx].reach_dist >= start_reach ) {
				end_idx--;
			}
			return std::make_pair( sda.begin_idx, end_idx + 1 );
		}
		assert( false );
		return { 0,0 };

	};
	const auto valid_combination = [&reach_dists, &chi, &min_pts, &get_reach_dist]( const SDA& sda, std::size_t sua_begin_idx, std::size_t sua_end_idx ) -> bool {
		if ( sda.mib > get_reach_dist( sua_end_idx + 1 ) * (1 - chi) ) { return false; }
		if ( sua_begin_idx - sda.end_idx < min_pts - 3 ) { return false; }
		//TODO: Checked conditions 1, 2, 3a?
		return true;
	};

	for ( std::size_t idx = 0; idx < reach_dists.size(); idx++ ) {
		double reach_i = reach_dists[idx].reach_dist;

		//Start of Steep Down Area?
		if ( idx < reach_dists.size() && is_steep_down_pt( idx ) ) {
			if ( reach_i > mib ) { mib = reach_i; }
			filter_sdas();
			std::size_t sda_end_idx = get_sda_end( idx );
			SDAs.push_back( SDA( idx, sda_end_idx, 0.0 ) );
			idx = sda_end_idx;
			if ( idx < reach_dists.size() - 1 ) { mib = reach_dists[idx + 1].reach_dist; }
			continue;
		}
		//Start of Steep Up Area?
		else if ( idx < reach_dists.size() && is_steep_up_pt( idx ) ) {
			filter_sdas();
			std::size_t sua_end_idx = get_sua_end( idx );

			for ( auto& sda : SDAs ) {
				if ( valid_combination( sda, idx, sua_end_idx ) ) {
					clusters.push_back( cluster_borders( sda, idx, sua_end_idx ) );
				}
			}
			idx = sua_end_idx;
			if ( idx < reach_dists.size() - 1 ) { mib = reach_dists[idx + 1].reach_dist; }
		}
		else { if ( reach_i > mib ) { mib = reach_i; } }
	}
	return clusters;
}


namespace internal{

std::vector<cluster_tree> flat_clusters_to_tree( const std::vector<chi_cluster_indices>& clusters_flat ){
	//sort clusters_flat such that children are ordered before their parents in clusters_flat_sorted
	std::vector<fplus::maybe<chi_cluster_indices>> clusters_flat_sorted_m( clusters_flat.size(), fplus::nothing<chi_cluster_indices>() );
	std::size_t next_free_idx = 0;
	for ( std::size_t idx = 0; idx < clusters_flat.size(); idx++ ) {
		while ( next_free_idx < clusters_flat_sorted_m.size() && clusters_flat_sorted_m[next_free_idx].is_just() ) {
			next_free_idx++;
		}
		std::size_t idx_pos = next_free_idx;
		std::size_t following_idx = idx + 1;
		while ( following_idx < clusters_flat.size() && clusters_flat[following_idx].second <= clusters_flat[idx].second ) {
			following_idx++;
			idx_pos++;
		}
		clusters_flat_sorted_m[idx_pos] = fplus::just( clusters_flat[idx] );
	}

	auto clusters_flat_sorted = fplus::justs( clusters_flat_sorted_m );
	assert( clusters_flat_sorted.size() == clusters_flat.size() );

	//compute tree from clusters_flat_sorted
	std::vector<cluster_tree> result;
	std::vector<cluster_tree> cluster_trees = fplus::transform( []( const chi_cluster_indices& c ) -> cluster_tree {
		return cluster_tree( c );
	}, clusters_flat_sorted );

	auto get_first_parent_idx = [&cluster_trees]( std::size_t idx )->std::size_t {
		auto cluster = cluster_trees[idx].get_root().get_data();
		for ( std::size_t first_parent_idx = idx + 1; first_parent_idx < cluster_trees.size(); first_parent_idx++ ) {
			auto parent_cluster = cluster_trees[first_parent_idx].get_root().get_data();
			if ( cluster.first >= parent_cluster.first && cluster.second <= parent_cluster.second ) {
				return first_parent_idx;
			}
		}
		return cluster_trees.size();
	};
	for ( std::size_t idx = 0; idx < cluster_trees.size(); idx++ ) {
		auto first_parent_idx = get_first_parent_idx( idx );
		if ( first_parent_idx >= cluster_trees.size() ) {
			result.push_back( cluster_trees[idx] );
		}
		else {
			cluster_trees[first_parent_idx].get_root().add_child( cluster_trees[idx].get_root() );
		}
	}

	return result;
}


inline void draw_cluster( bgr_image& cluster_indicator_img, const Node<chi_cluster_indices>& cluster, const std::size_t depth,
						  const double x_norm, const std::size_t v_dist ) {
	std::size_t v_offset = (depth+1)*v_dist;
	assert( v_offset <= cluster_indicator_img.size().height_ );
	std::size_t height = cluster_indicator_img.size().height_ - v_offset;
	img_pos x1( fplus::round<double, std::size_t>(x_norm*cluster.get_data().first), height );
	img_pos x2( fplus::round<double, std::size_t>( x_norm*cluster.get_data().second), height );
	bgr_col col( 0, 0, 0 );
	bgr_col st_col( 0, 255, 0 );
	bgr_col end_col( 255, 0, 0 );
	plot_line_segment( cluster_indicator_img, x1, x2, col );
	plot_pixel( cluster_indicator_img, x1, st_col );
	plot_pixel( cluster_indicator_img, x2, end_col );
	cluster_indicator_img.save( "./tmp" );
	for ( const auto& c : cluster.get_children() ) {
		draw_cluster( cluster_indicator_img, c, depth+1, x_norm, v_dist );
	}
}

}//namepsace internal


inline std::vector<cluster_tree> get_chi_clusters( const std::vector<reachability_dist>& reach_dists, const double chi, std::size_t min_pts ) {
	auto clusters_flat = get_chi_clusters_flat( reach_dists, chi, min_pts );
	return internal::flat_clusters_to_tree( clusters_flat );
}


inline bgr_image  draw_reachability_plot_with_chi_clusters( const std::vector<reachability_dist>& reach_dists,
															const double chi, const std::size_t min_pts,
															const std::size_t min_width = 100)
{
	auto img = draw_reachability_plot( reach_dists, min_width );
	auto cluster_trees = optics::get_chi_clusters( reach_dists, chi, min_pts );

	std::size_t max_tree_depth = 0;
	for ( const auto& t : cluster_trees ) {
		auto depth = optics::tree_depth<chi_cluster_indices>( t.get_root() );
		if ( depth > max_tree_depth ) { max_tree_depth = depth; }
	}
	std::size_t v_space = 2;
	bgr_image cluster_indicator_img( size_2d( img.size().width_, (max_tree_depth +1) * v_space ), bgr_col(255,255,255) );
	
	double x_norm = 1;
	if ( min_width > reach_dists.size() ) {
		x_norm =  static_cast<double>(min_width)/ static_cast<double>(reach_dists.size()-1);
	}
	for ( const auto& t : cluster_trees ) {
		internal::draw_cluster( cluster_indicator_img, t.get_root(), 0, x_norm, v_space );
	}
	cluster_indicator_img.save( "./cluster_indicator" );
	img.append_rows(cluster_indicator_img);
	return img;
}


template<typename T>
bgr_image draw_2d_clusters( const std::vector<std::vector<geom::Vec<T,2>>>& clusters ) {
	auto box = geom2d::bounding_box( fplus::concat( clusters ) );
	bgr_image cluster_image( size_2d( fplus::round<double, std::size_t>(box.get_size().first+1), fplus::round<double, std::size_t>( box.get_size().second+1)), bgr_col(255,255,255) );
	std::array<bgr_col, 6> colours = { bgr_col( 255,0,0 ), bgr_col( 0,255,0 ), bgr_col(0,0,255), bgr_col(255,255,0), bgr_col(255,0,255), bgr_col(0,255,255) };
	int col_idx = 0;
	for ( const auto& cluster : clusters ) {
		bgr_col col = colours[col_idx];
		++col_idx %= colours.size();
		auto cluster_box = geom2d::bounding_box( cluster );
		for( const auto& edge: fplus::overlapping_pairs_cyclic(cluster_box.points()) ){
            plot_line_segment( cluster_image,
							   img_pos( fplus::round<double, std::size_t>(edge.first.x() - box.bl().x() ), fplus::round<double, std::size_t>( edge.first.y() - box.bl().y() )),
							   img_pos( fplus::round<double, std::size_t>( edge.second.x() - box.bl().x() ), fplus::round<double, std::size_t>( edge.second.y() - box.bl().y() ) ),
							   col );
		}
		for ( const auto& pt : cluster ) {
			img_pos cluster_pt = img_pos( fplus::round<double, std::size_t>( pt.x() - box.bl().x() ), fplus::round<double, std::size_t>( pt.y() - box.bl().y() ) );
			plot_circle( cluster_image, cluster_pt, 2, col );
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

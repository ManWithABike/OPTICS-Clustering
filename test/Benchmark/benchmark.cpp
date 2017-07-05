// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)


#include "../../include/optics/optics.hpp"
#include "../../include/optics/Stopwatch.hpp"

#include <random>

namespace sw = stopwatch;


//static const std::random_device rd;
static std::mt19937 gen( 1 );


template <typename T, typename = typename std::enable_if<std::is_integral<T>::value >::type >
std::uniform_int_distribution<T> uniform_distribution( T a, T b )
{
	return std::uniform_int_distribution<T>( a, b );
}

template <typename T, typename = typename std::enable_if<std::is_floating_point<T>::value >::type >
std::uniform_real_distribution<T> uniform_distribution( T a, T b )
{
	return std::uniform_real_distribution<T>( a, b );
}


template<std::size_t dimension, typename type, typename distribution>
std::array<type, dimension> get_random_point( distribution& dis ) {
	std::array<type, dimension> result;
	for ( std::size_t d = 0; d < dimension; d++ ) {
		result[d] = dis( gen );
	}
	return result;
}


//Returns the estimated space for each point given a perfect uniform distribution of n_points in dim dimensional volume space_volume
double eps_guess(double space_volume, std::size_t n_points, std::size_t dim ) {
	double d = static_cast<double> (dim);
	double nominator = space_volume * 1.0 * std::tgamma( d / 2.0 + 1.0 );
	double denominator = static_cast<double>(n_points) * std::sqrt( std::pow( geom::pi, d ) );
	double r = std::pow( nominator / denominator, 1.0 / d );
	return r;
}


template<std::size_t n_points, std::size_t dimension, std::uint64_t space_volume, typename type, std::size_t laps>
std::uint64_t test( std::size_t min_pts, double epsilon = -1.0 ) {
	static_assert(dimension > 0, "Dimension must be >0");

	//Space parameters
	double edge_length = std::pow( static_cast<double>(space_volume), 1.0 / static_cast<double>(dimension) );
	std::cout << "Space volume: " << std::pow(edge_length, static_cast<double>(dimension) ) << std::endl;
	std::cout << "Space edge length: " << edge_length << std::endl;
	double eps = eps_guess( static_cast<double>(space_volume), n_points, dimension );
	std::cout << "Estimated space per point: " << eps << std::endl;
	
	//Create a distribution that outputs coordinates between 0 and edge_length

	//std::uniform_real_distribution<type> dis( 0, edge_length );
	auto dis = uniform_distribution<type>( 0, static_cast<type>(edge_length) );

	//Create n_points random points
	std::vector<std::array<type, dimension>> points;
	points.reserve( n_points );
	for ( std::size_t n = 1; n <= n_points; n++ ) {
		points.push_back( get_random_point<dimension, type>( dis ) );
	}

	//Let's go
	std::cout << "Start " << laps << " computations of  optics::compute_reachability_dist()" << std::endl;
	sw::Stopwatch watch;
	for ( std::size_t lap = 1; lap <= laps; lap++ ) {
		if ( lap % (laps/10) == 0 ) std::cout << lap << "..";
		optics::compute_reachability_dists( points, min_pts, epsilon );
		watch.lap();
	}

	auto elapsed = watch.elapsed_laps();
	std::cout << std::endl;
	std::cout << "Lap times [ms]: " << sw::show_times( elapsed.second ) << std::endl;
	std::cout << "Total time: " << elapsed.first << "ms" << std::endl;
	std::cout << "Mean time: " << elapsed.first/laps << "ms" << std::endl;

	return elapsed.first / laps;
}


int main() {
	std::cout << "OPTICS Benchmark" << std::endl;

	test<100000, 2, 1000*1000, double, 10>( 10 );

	return 0;
}
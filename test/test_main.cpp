#include "../include/optics/optics.h"
#include <vector>



void test_1(){
	static const int N = 2;
	typedef std::array<double, N> point; //A list of N cartesian coordinates makes a point

    std::vector<point> points; //Your list of points goes here
	points = { {100,100}, {102,100}, {101,101},           //cluster 1
			   {-1,0}, {1,0}, {0,1},                     //cluster 2
			   {-100,-100}, {-102,-100}, {-101,-101}     //cluster 3
	};
	auto reach_dists = optics::compute_reachability_dists( points, 2, 10 );
	for( const auto& x : reach_dists){
        std::cout << x.to_string() << "; ";
	}

	auto clusters = optics::get_cluster_indices(reach_dists, 10);
	assert(clusters.size() == 3);
	assert( ( fplus::sort( clusters[0] ) == std::vector <std::size_t>({ 0,1,2 }) ) );
	assert( ( fplus::sort( clusters[1] ) == std::vector <std::size_t>({ 3,4,5 }) ) );
	assert( ( fplus::sort( clusters[2] ) == std::vector <std::size_t>({ 6,7,8 }) ) );

	optics::draw_reachability_plot( reach_dists, "ReachabilityPlot" );
}

void test_2() {
	static const int N = 2;
	typedef std::array<int, N> point; //A list of N cartesian coordinates makes a point

	std::vector<point> points; //Your list of points goes here
	points = { { 100,100 },{ 102,100 },{ 101,101 },           //cluster 1
	{ -1,0 },{ 1,0 },{ 0,1 },                     //cluster 2
	{ -100,-100 },{ -102,-100 },{ -101,-101 }     //cluster 3
	};
	auto reach_dists = optics::compute_reachability_dists( points, 2, 10 );
	for ( const auto& x : reach_dists ) {
		std::cout << x.to_string() << "; ";
	}

	auto clusters = optics::get_cluster_indices( reach_dists, 10 );
	assert( clusters.size() == 3 );
	assert( (fplus::sort( clusters[0] ) == std::vector <std::size_t>( { 0,1,2 } )) );
	assert( (fplus::sort( clusters[1] ) == std::vector <std::size_t>( { 3,4,5 } )) );
	assert( (fplus::sort( clusters[2] ) == std::vector <std::size_t>( { 6,7,8 } )) );

	optics::draw_reachability_plot( reach_dists, "ReachabilityPlot" );
}


int main()
{
    test_1();

	test_2();
}

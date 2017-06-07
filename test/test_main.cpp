#include <optics/optics.h>
#include <vector>

static const int N = 2;
typedef std::array<double, N> point; //A list of N cartesian coordinates makes a point

int main()
{
	std::vector<point> points; //Your list of points goes here
	points = { {100,100}, {102,100}, {101,101},           //cluster 1
			   {-1,0}, {1,0}, {0,1},                     //cluster 2
			   {-100,-100}, {-102,-100}, {-101,-101}     //cluster 3
	};
	auto reach_dists = optics::compute_reachability_dists<double,2>( points, 2, 10 );
	for( const auto& x : reach_dists){
        std::cout << x.to_string() << "; ";
	}
	optics::draw_reachability_plot( reach_dists, "reachdists.bmp" );

}

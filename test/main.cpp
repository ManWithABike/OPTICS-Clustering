#include <optics.h>

typedef std::vector<double> point; //A list of n cartesian coordinates makes a point
std::vector<point> points; //Your list of points goes here

int main(){
   auto reach_dists = optics::compute_reachability_dists( points, 10, 100 );
   optics::draw_reachability_plot( reach_dists, "D:/reachdists.bmp" );
}

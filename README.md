# OPTICS-Clustering

**Ordering points to identify the clustering structure (OPTICS)** is an algorithm for finding density-based clusters in spatial data. It was [presented](https://github.com/CrikeeIP/OPTICS-Clustering/blob/master/background/OPTICS.pdf) by Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel and Jörg Sander in 1999.

##Introduction
This repository is home to a **C++** implementation of the OPTICS algorithm as described by Ankherst et al. .
For further explanation on how the algorithm works, see e.g. [Wikipedia](https://en.wikipedia.org/wiki/OPTICS_algorithm) or [YouTube}(https://www.youtube.com/watch?v=8kJjgowewOs).
It relies on the boost RTree implementation in order to efficiently find neighbourhoods of a given point.


##Usage
Suppose you have a set of points in R^n, described in cartesion coordinates, and wonder if they have a cluster structure.
Then you might consider using this library, as it offers an interface that lets you draw a [reachability-plot]() with two lines of code:

```cpp
#include<optics.h>

typedef std::vector<double> point; //A list of n cartesian coordinates makes a point
std::vector<point> points; //Your list of points goes here

int main(){
   auto reach_dists = optics::compute_reachability_dists( points, 10, 100 );
   optics::draw_reachability_plot( reach_dists, "D:/reachdists.bmp" );
}
```


##Dependencies
Three lightweight header-only libraries:  
1. [Geometry](https://github.com/CrikeeIP/Geometry)  
2. [FunctionalPlus](https://github.com/Dobiasd/FunctionalPlus)  
3. [CImg](https://github.com/dtschump/CImg)

And boost (::geometry and ::index, to be exact)  
4. [Boost](http://www.boost.org/)

In order to do so, the points of the database are (linearly) ordered such that points which are spatially closest become neighbors in the ordering. Additionally, a special distance is stored for each point that represents the density that needs to be accepted for a cluster in order to have both points belong to the same cluster.

##Disclaimer

The functionality in this library initially grew due to my personal need for it while using C++ on a regular basis. I try my best to make it error free and as comfortable to use as I can. The API still might change in the future. If you have any suggestions, find errors, miss some functions or want to give general feedback/criticism, I'd love to hear from you. Of course, contributions are also very welcome.

##License

Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

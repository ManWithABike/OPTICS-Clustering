[![Build Status](https://travis-ci.org/CrikeeIP/OPTICS-Clustering.svg?branch=master)][travis]

[travis]: https://travis-ci.org/CrikeeIP/OPTICS-Clustering

# OPTICS-Clustering (UNDER CONSTRUCTION)

**Ordering points to identify the clustering structure ([OPTICS](https://github.com/CrikeeIP/OPTICS-Clustering/blob/master/background/OPTICS.pdf))** is an algorithm for finding density-based clusters in spatial data. It was presented by Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel and Jörg Sander in 1999.

## Introduction
This repository is home to a C++ implementation of the OPTICS algorithm as described by Ankerst et al. . It aims at providing an easy-to-use clustering algorithm which does not require knowledge of the number of clusters a priori.  
For further explanation on how the algorithm works, see e.g. [Wikipedia](https://en.wikipedia.org/wiki/OPTICS_algorithm) or [YouTube](https://www.youtube.com/watch?v=8kJjgowewOs).  
The implementation relies on the `boost::rtree` to efficiently retrieve the neighbourhood of a given point.


## Usage
Suppose you have a set of points in R^n, described in Cartesian coordinates, and wonder if they have a cluster structure.
Then you might consider using this library, as it offers an interface that lets you visually inspect the cluster structure of the data space (using a [reachability-plot](https://github.com/CrikeeIP/OPTICS-Clustering/blob/master/resources/reachabilityplot.png)) - *and* extract those clusters with three lines of code:

```cpp
#include <optics/optics.h>

typedef std::vector<double> point; //A list of n cartesian coordinates makes a point
std::vector<point> points; //Your list of points goes here

int main(){
   auto reach_dists = optics::compute_reachability_dists<T,N>( points, min_pts, epsilon );
   optics::draw_reachability_plot( reach_dists, "./reachdists.bmp" );
   auto clusters = optics::get_cluster_indices( reach_dists, threshold );
}
```


## Dependencies
Two lightweight header-only libraries, and boost (`::geometry` and `::index`, to be exact) 
1. [Geometry](https://github.com/CrikeeIP/Geometry)  
2. [FunctionalPlus](https://github.com/Dobiasd/FunctionalPlus)  
3. [Boost](http://www.boost.org/)


## Installation
Before the installation, make sure you have installed [Boost](http://www.boost.org/).
Subsequently, you can use one of the following two alternative ways to install OPTICS-Clustering:

***Alternative 1:***  
Clone this repository
```sh
git clone https://github.com/CrikeeIP/OPTICS-Clustering
```
and execute the the `install.sh` script delivered with it.

***Alternative 2:***  
If you're uncomfortable running who-knows-what-they'll-do foreign scripts (and are too tired to check them before execution), you can also do it manually:
1. Clone & include the [Geometry](https://github.com/CrikeeIP/Geometry) header-lib,
2. Clone & install [FunctionalPlus](https://github.com/Dobiasd/FunctionalPlus) 
3. Clone & install the OPTICS repository:
```sh
git clone https://github.com/CrikeeIP/OPTICS-Clustering
cd OPTICS-clustering
./configure && make && make install
```

## Disclaimer

This librarys functionality initially grew due to my personal need for it - an easy to use clustering algorithm which does not need to know the number of clusters a priori.
I try my best to make it error free and as comfortable to use as I can. The API still might change in the future. If you have any suggestions, find errors, miss some functions or want to give general feedback/criticism, I'd love to hear from you. Of course, [contributions](https://github.com/CrikeeIP/OPTICS-Clustering/pulls) are also very welcome.

## License

Distributed under the MIT Software License (X11 license). (See accompanying file LICENSE.)

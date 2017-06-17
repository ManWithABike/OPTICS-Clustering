[![Build Status](https://travis-ci.org/CrikeeIP/OPTICS-Clustering.svg?branch=master)][travis]

[travis]: https://travis-ci.org/CrikeeIP/OPTICS-Clustering

# OPTICS-Clustering (UNDER CONSTRUCTION)

**Ordering points to identify the clustering structure (OPTICS)** is an algorithm for finding density-based clusters in spatial data. It was [presented](https://github.com/CrikeeIP/OPTICS-Clustering/blob/master/background/OPTICS.pdf) by Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel and Jörg Sander in 1999.

## Introduction
This repository is home to a C++ implementation of the OPTICS algorithm as described by Ankherst et al. . It aims at providing an easy-to-use clustering algorithm which does not require knowledge of the number of clusters a priori.  
For further explanation on how the algorithm works, see e.g. [Wikipedia](https://en.wikipedia.org/wiki/OPTICS_algorithm) or [YouTube](https://www.youtube.com/watch?v=8kJjgowewOs).
The implementation relies on the Boost RTree for efficient neighbourhood queries of a given point.


## Usage
Suppose you have a set of points in R^n, described in cartesian coordinates, and wonder if they have a cluster structure.
Then you might consider using this library, as it offers an interface that lets you extract clusters *and* draw the corresponding [reachability-plot](https://github.com/CrikeeIP/OPTICS-Clustering/blob/master/resources/reachabilityplot.png) with three lines of code:

```cpp
#include <optics/optics.h>

typedef std::vector<double> point; //A list of n cartesian coordinates makes a point
std::vector<point> points; //Your list of points goes here

int main(){
   auto reach_dists = optics::compute_reachability_dists<T,N>( points, min_pts, epsilon );
   optics::make_clusters( reach_dists, threshold );
   optics::draw_reachability_plot( reach_dists, "./reachdists.bmp" );
}
```


## Dependencies
Two lightweight header-only libraries:  
1. [Geometry](https://github.com/CrikeeIP/Geometry)  
2. [FunctionalPlus](https://github.com/Dobiasd/FunctionalPlus)  
   
   And boost (`::geometry` and `::index`, to be exact) 
3. [Boost](http://www.boost.org/)


## Installation
Before the installation, make sure you have installed [Boost](http://www.boost.org/).
Subsequently, you can use one of the following two alternative ways to install OPTICS-Clustering:

***Alternative 1:***  
Clone this repository
```sh
git clone https://github.com/CrikeeIP/OPTICS-Clustering
```
and execute the the *install.sh* script delivered with it.

***Alternative 2:***  
If you're uncomfortable running who-knows-what-they'll-do foreign scripts (and are too tired to check them before execution), you can also do it manually:
```sh
#Clone the reporitory
git clone https://github.com/CrikeeIP/OPTICS-Clustering
#Download CImg header
cd OPTICS-Clustering
cd include
cd optics
mkdir CImg
cd Cimg
wget https://raw.githubusercontent.com/dtschump/CImg/master/CImg.h
cd ..
#Download Geometry header
mkdir Geometry
cd Geometry
wget https://raw.githubusercontent.com/CrikeeIP/Geometry/master/include/geometry/geometry.h
cd ..
cd ..
cd ..

#Clone the FunctionalPlus repository and install it
git clone https://github.com/Dobiasd/FunctionalPlus
cd FunctionalPlus
mkdir build
cd build
cmake ..
sudo make install
cd ..
cd ..

#Run the test
cd OPTICS-Clustering
cd test
g++ --std=c++11 -I../include main.cpp -lX11 -lpthread
./a.out
```


## Disclaimer

This librarys functionality initially grew due to my personal need for it - an easy to use clustering algorithm which does not need to know the number of clusters a priori.
I try my best to make it error free and as comfortable to use as I can. The API still might change in the future. If you have any suggestions, find errors, miss some functions or want to give general feedback/criticism, I'd love to hear from you. Of course, [contributions](https://github.com/CrikeeIP/OPTICS-Clustering/pulls) are also very welcome.

## License

Distributed under the MIT Software License (X11 license). (See accompanying file LICENSE.)

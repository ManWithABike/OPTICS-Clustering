# OPTICS-Clustering

**Ordering points to identify the clustering structure (OPTICS)** is an algorithm for finding density-based clusters in spatial data. It was [presented](https://github.com/CrikeeIP/OPTICS-Clustering/blob/master/background/OPTICS.pdf) by Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel and Jörg Sander in 1999.

##Introduction
This repository is home to a **C++** implementation of the OPTICS algorithm as described by Ankherst et al. .
For further explanation on how the algorithm works, see e.g. [Wikipedia](https://en.wikipedia.org/wiki/OPTICS_algorithm) or [YouTube](https://www.youtube.com/watch?v=8kJjgowewOs).
It relies on the boost RTree implementation in order to efficiently find neighbourhoods of a given point.


##Usage
Suppose you have a set of points in R^n, described in cartesion coordinates, and wonder if they have a cluster structure.
Then you might consider using this library, as it offers an interface that lets you draw a [reachability-plot](https://github.com/CrikeeIP/OPTICS-Clustering/blob/master/resources/reachabilityplot.png) with two lines of code:

```cpp
#include <optics.h>

typedef std::vector<double> point; //A list of n cartesian coordinates makes a point
std::vector<point> points; //Your list of points goes here

int main(){
   auto reach_dists = optics::compute_reachability_dists( points, 5, 40 );
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


##Installation
Before the installation, make sure you have installed Boost.
Subsequently, you can use one of the following two alternative ways to install OPTICS-Clustering:

***Alternative 1:***
Clone this repository
```sh
git clone https://github.com/CrikeeIP/OPTICS-Clustering
```
and execute the the *install.sh* script delivered with it.

***Alternative 2:***
If you're uncomfortable running who-knows-what-they'll-do foreign scripts (and are too tired to check them before executing), you can do it manually:
```sh
#Clone the reporitory
git clone https://github.com/CrikeeIP/OPTICS-Clustering
#Download CImg header
cd OPTICS-Clustering
cd include
cd optics
mkdir CImg
wget https://raw.githubusercontent.com/dtschump/CImg/master/CImg.h
cd ..
#Download Geometry header
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
gcc -I../include main.cpp
./a.out
```


##Disclaimer

The functionality in this library initially grew due to my personal need for it while using C++ on a regular basis. I try my best to make it error free and as comfortable to use as I can. The API still might change in the future. If you have any suggestions, find errors, miss some functions or want to give general feedback/criticism, I'd love to hear from you. Of course, [contributions](https://github.com/CrikeeIP/OPTICS-Clustering/pulls) are also very welcome.

##License

Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

# OPTICS-Clustering

**Ordering points to identify the clustering structure (OPTICS)** is an algorithm for finding density-based clusters in spatial data. It was [presented](https://github.com/CrikeeIP/OPTICS-Clustering/blob/master/background/OPTICS.pdf) by Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel and Jörg Sander in 1999.

##Introduction
This repository is home to a **C++** implementation of the OPTICS algorithm as described by Ankherst et al. .


##Dependencies
[FunctionalPlus](https://github.com/Dobiasd/FunctionalPlus)
[CImg](https://github.com/dtschump/CImg)
[Boost](http://www.boost.org/)

In order to do so, the points of the database are (linearly) ordered such that points which are spatially closest become neighbors in the ordering. Additionally, a special distance is stored for each point that represents the density that needs to be accepted for a cluster in order to have both points belong to the same cluster.

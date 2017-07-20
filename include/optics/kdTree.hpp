// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)


#pragma once

#include <array>
#include <fplus/fplus.hpp>

template<typename T, std::size_t dimension, std::size_t n_points>
class KDTree {
public:
	typedef std::array<T, dimension> point_type;

	KDTree( const std::array<point_type, n_points>& points ) : points_(points)
	{
		std::array<std::size_t, n_points> indices = fplus::numbers<std::array<std::size_t, n_points>>(0, n_points);
		
		std::_Nth_element_unchecked()

	}


private:
	point_type& points_;
};
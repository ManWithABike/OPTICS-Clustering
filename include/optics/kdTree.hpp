// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)


#pragma once

#include <array>
#include <algorithm>
#include <fplus/fplus.hpp>



namespace kdt{

template<std::size_t N>
using Vec = std::array<double, N>;

template<std::size_t N>
double sum( const Vec<N>& x ) {
	size_t i = 0;
	//TODO: AVX2 loop

	/*
	//TODO: replacement for _mm256_cvtsd_f64
	//AVX loop
	__m256d vsum0_avx = _mm256_setzero_pd();
	__m256d vsum1_avx = _mm256_setzero_pd();
	for ( ; i < (size & ~0x3); i += 4 )
	{
	__m256d v0 = _mm256_load_pd( &x[i] );
	__m256d v1 = _mm256_load_pd( &x[i + 4] );
	vsum0_avx = _mm256_add_pd( vsum0_avx, v0 );
	vsum1_avx = _mm256_add_pd( vsum1_avx, v1 );
	}
	vsum0_avx = _mm256_add_pd( vsum0_avx, vsum1_avx );    // vertical ops down to one accumulator
	vsum0_avx = _mm256_hadd_pd( vsum0_avx, vsum0_avx );   // horizontal add of the single register
	double sum1 = _mm256_cvtsd_f64( vsum0_avx ); //TODO: replacement for _mm256_cvtsd_f64
	*/

	//SSE2 loop
	__m128d vsum0 = _mm_setzero_pd();
	__m128d vsum1 = _mm_setzero_pd();
	for ( ; i < (N & ~0x3); i += 4 )
	{
		__m128d v0 = _mm_load_pd( &x[i] );
		__m128d v1 = _mm_load_pd( &x[i + 2] );
		vsum0 = _mm_add_pd( vsum0, v0 );
		vsum1 = _mm_add_pd( vsum1, v1 );
	}
	vsum0 = _mm_add_pd( vsum0, vsum1 );    // vertical ops down to one accumulator
	vsum0 = _mm_hadd_pd( vsum0, vsum0 );   // horizontal add of the single register
	double sum2 = _mm_cvtsd_f64( vsum0 );

	double sum = sum2; //+sum1;
						// Serial loop
	for ( ; i < N; i++ )
	{
		sum += x[i];
	}
	return sum;
}

template<std::size_t N>
Vec<N> subtract( const Vec<N>& x, const Vec<N>& y ) {
	//void add( double* result, const double* a, const double* b, size_t size )
	size_t i = 0;
	Vec<N> result;
	// Note we are doing as many blocks of 8 as we can.  If the size is not divisible by 8
	// then we will have some left over that will then be performed serially.
	// AVX-512 loop
	/*for ( ; i < (size & ~0x7); i += 8 )
	{
	const __m512d kA8 = _mm512_load_pd( &x[i] );
	const __m512d kB8 = _mm512_load_pd( &y[i] );

	const __m512d kRes = _mm512_add_pd( kA8, kB8 );
	_mm512_store_pd( &result[i], kRes );
	}*/

	auto x_ptr = &x[0];
	// AVX loop
	for ( ; i < (N & ~0x3); i += 4 )
	{
		
		const __m256d kA4 = _mm256_load_pd( &x[i] );
		const __m256d kB4 = _mm256_load_pd( &y[i] );

		const __m256d kRes = _mm256_sub_pd( kA4, kB4 );
		_mm256_store_pd( &result[i], kRes );
	}

	// SSE2 loop
	for ( ; i < (N & ~0x1); i += 2 )
	{
		const __m128d kA2 = _mm_load_pd( &x[i] );
		const __m128d kB2 = _mm_load_pd( &y[i] );

		const __m128d kRes = _mm_sub_pd( kA2, kB2 );
		_mm_store_pd( &result[i], kRes );
	}

	// Serial loop
	for ( ; i < N; i++ )
	{
		result[i] = x[i] - y[i];
	}
	return result;
}


template<std::size_t N>
Vec<N> product( const Vec<N>& x, const Vec<N>& y ) {
	size_t i( 0 );
	Vec<N> result;

	// Note we are doing as many blocks of 8 as we can.  If the size is not divisible by 8
	// then we will have some left over that will then be performed serially.
	// AVX-512 loop
	/*for ( ; i < (size & ~0x7); i += 8 )
	{
	const __m512d kA8 = _mm512_load_pd( &x[i] );
	const __m512d kB8 = _mm512_load_pd( &y[i] );

	const __m512d kRes = _mm512_mul_pd( kA8, kB8 );
	_mm512_store_pd( &result[i], kRes );
	}*/

	for ( ; i < (N & ~0x3); i += 4 )
	{
		const __m256d kX4 = _mm256_load_pd( &x[i] );
		const __m256d kY4 = _mm256_load_pd( &y[i] );

		const __m256d kRes = _mm256_mul_pd( kX4, kY4 );
		_mm256_store_pd( &result[i], kRes );
	}

	// SSE2 loop
	for ( ; i < (N & ~0x1); i += 2 )
	{
		const __m128d kX2 = _mm_load_pd( &x[i] );
		const __m128d kY2 = _mm_load_pd( &y[i] );

		const __m128d kRes = _mm_mul_pd( kX2, kY2 );
		_mm_store_pd( &result[i], kRes );
	}

	// Serial loop
	for ( ; i < N; i++ )
	{
		result[i] = x[i] * y[i];
	}
	return result;
}


template<std::size_t N>
__declspec(noinline) double square_distance_intrinsics( const Vec<N>& x, const Vec<N>& y ) {
	auto diff = subtract( x, y );
	return sum( product( diff, diff ) );
}


constexpr bool is_powerof2(std::size_t v) {
	return v && ((v & (v - 1)) == 0);
}

enum DistanceTypes { MANHATTAN, EUCLIDEAN, MAXIMUM };

template<typename CoordsType, std::size_t dimension, std::size_t n_points_total, std::size_t n_points, std::size_t max_points_per_node, std::size_t split_dim, typename Parent>class KDTree;

template<typename CoordsType, std::size_t dimension>
static double square_distance( const std::array<CoordsType, dimension>& p1, const std::array<CoordsType, dimension>& p2 ) {
	double result = 0.0;
	//std::array<CoordsType, dimension> temp;
	for ( std::size_t i = 0; i < dimension; i++ ) {
		//temp[i] = (p1[i] - p2[i]);
		double d = (p1[i] - p2[i]);
		result +=  d*d;
	}
	//for ( std::size_t i = 0; i < dimension; i++ ) {
//		result += temp[i];
//	}
	return result;
}


template<typename CoordsType, std::size_t dimension, std::size_t n_points_total, std::size_t n_points, typename Parent>
class KDTreeLeaf {
	friend Parent;
	/*
	template<CoordsType, dimension, 2 * n_point, std::size_t max_points_per_node, std::size_t split_dim>
	friend class KDTree;
	template<CoordsType, dimension, 2 * n_points + 1, std::size_t max_points_per_node, std::size_t split_dim>
	friend class KDTree;
	*/
public:
	typedef std::array<CoordsType, dimension> Point;
	KDTreeLeaf() {}
	KDTreeLeaf( const std::vector<Point>& points_,
				const std::array<std::size_t, n_points>& indices_
	) 
		: indices( indices_ ), points( &points_ )
	{
		if ( points_.size() != n_points_total ) {
			std::cerr << "KDTreeLeaf(): points_ does not contain n_points_total elements! Abort";
			std::exit( 1 );
		}
	}

	void radius_search_ ( const Point& p, double radius, std::vector<std::size_t>& neighbors ) const //TODO: make private & friend KDTree
	{
		for ( const auto& i : indices ) {
			if ( square_distance_intrinsics( (*points)[i], p ) <= radius*radius ) {
				neighbors.push_back( i );
			}
		}
	}

private:
	std::array<std::size_t, n_points> indices;
	std::vector<Point> const* points;
};


/*
template<typename CoordsType, std::size_t dimension, std::size_t n_points, std::size_t max_points_per_node, std::size_t split_dim>
std::enable_if_t<
	(n_points <= max_points_per_node),
	std::unique_ptr<KDTree<CoordsType, dimension, n_points, max_points_per_node, (split_dim + 1) % dimension>>
>
make_child( const std::array<std::array<CoordsType, dimension>, n_points>& points_ ) {
	return std::make_unique<KDTreeLeaf<CoordsType, dimension, n_points, max_points_per_node, (split_dim + 1) % dimension>>( points_ );
}

template<typename CoordsType, std::size_t dimension, std::size_t n_points, std::size_t max_points_per_node, std::size_t split_dim>
std::enable_if_t<
	(n_points > max_points_per_node),
	std::unique_ptr<KDTree<CoordsType, dimension, n_points, max_points_per_node, (split_dim + 1) % dimension>>
>
make_child( const std::array<std::array<CoordsType, dimension>, n_points>& points_ ) {
	return std::make_unique<KDTree<CoordsType, dimension, n_points, max_points_per_node, (split_dim + 1) % dimension>>( points_ );
}
*/

template <typename T, std::size_t N>
std::unique_ptr<std::array<T, N>> make_array()
{
	return std::make_unique<std::array<T, N>>();
}


template<std::size_t N>
constexpr std::unique_ptr<std::array<std::size_t, N>> indices() {
	auto result = make_array<std::size_t, N>();
	for ( std::size_t i = 0; i < N; i++ ) {
		(*result)[i] = i;
	}
	return result;
}


template<typename CoordsType, std::size_t dimension, std::size_t n_points, std::size_t max_points_per_node>
std::unique_ptr<KDTree<CoordsType, dimension, n_points, n_points, max_points_per_node, 0, std::true_type>> 
make_KDTree( const std::vector<std::array<CoordsType, dimension>>& points_ )
{
	if ( points_.size() != n_points ) {
		std::cerr << "make_KDTree(): points_ does not contain n_points elements. Aborting.";
		std::exit( 1 );
	}
	//const auto idxs = indices<n_points>(); // todo unten rein
	//return std::make_unique<KDTree<CoordsType, dimension, n_points, n_points, max_points_per_node, 0, std::true_type>>( points_, idxs );
	//return std::unique_ptr<KDTree<CoordsType, dimension, n_points, n_points, max_points_per_node, 0, std::true_type>( new KDTree<CoordsType, dimension, n_points, n_points, max_points_per_node, 0, std::true_type>( points_, idxs ) );

	const auto idxs = indices<n_points>();
	std::unique_ptr<KDTree<CoordsType, dimension, n_points, n_points, max_points_per_node, 0, std::true_type>> res;
	auto naked_ptr = new KDTree<CoordsType, dimension, n_points, n_points, max_points_per_node, 0, std::true_type>( points_, *idxs );
	res.reset( naked_ptr );
	return res;
}


template<typename CoordsType, std::size_t dimension, std::size_t n_points_total, std::size_t n_points, std::size_t max_points_per_node, std::size_t split_dim, typename Parent>
class KDTree {
	//static_assert(is_powerof2(n_points), "The total number of points stored in a KDTree 'n_points' has to be a power of two.");
	//static_assert(is_powerof2(max_pts_per_node), "");

public:
	typedef std::array<CoordsType, dimension> Point;

	using Me = KDTree<CoordsType, dimension, n_points_total, n_points, max_points_per_node, split_dim, Parent>;

	friend std::unique_ptr<KDTree<CoordsType, dimension,
							n_points, n_points,
							max_points_per_node, 0, std::true_type>>
			make_KDTree<CoordsType, dimension, n_points, max_points_per_node>
			( const std::vector<std::array<CoordsType, dimension>>& points_ );
			
	//friend std::make_unique<Me, const std::array<std::array<CoordsType, dimension>, n_points>&, const std::array<std::size_t, n_points>&>( const std::array<std::array<CoordsType, dimension>, n_points>&, const std::array<std::size_t, n_points>& );
	
	std::vector<std::size_t> radius_search( const Point& p, double radius ) const
	{
		std::vector<std::size_t> neighbors;
		neighbors.reserve( 2 * max_points_per_node ); //TODO: intelligent guess of number of expected neighbors?
		radius_search_( p, radius, neighbors );
		neighbors.shrink_to_fit();
		return neighbors;
	}

	KDTree() {}

private:
	KDTree( const std::vector<Point>& points_, std::array<std::size_t, n_points>& indices_ )
	{
		if ( points_.size() != n_points_total ) {
			std::cerr << "KDTree(): points_ does not contain n_points_total elements! Aborting";
			std::exit( 1 );
		}
		//std::array<std::size_t, n_points> indices = fplus::numbers<std::size_t, std::array<std::size_t, n_points>>(static_cast<std::size_t>(0), n_points);
		//fplus::nth_element 
		/*auto split_index_it = std::begin(indices);
		std::advance(split_index_it, static_cast<typename split_index_it::difference_type>(n));
		auto split_index_it = std::nth_element(std::begin(indices), split_index_it, std::end(indices), comp<split_dim>);
		auto split_index = *split_index_it;
		auto split_point = points[split_index];
		split = split_point[split_dim];
		*/
		const auto comp = [&points_]( std::size_t i1, std::size_t i2 ) {
			return points_[i1][split_dim] < points_[i2][split_dim];
		};

		//auto split_index = fplus::nth_element_by(comp, n_points / 2, indices_ );
		auto middle = std::begin( indices_ );
		std::advance( middle, n_points / 2 );
		std::nth_element( std::begin( indices_ ), middle, std::end( indices_ ), comp );
		auto split_index = *middle;
		auto split_point = points_[split_index];
		split = split_point[split_dim];

		auto left_indices = make_array<std::size_t, n_points/2>(); //Could hold references instead of real points?
		auto right_indices = make_array<std::size_t, n_points - n_points / 2>();
		std::size_t l_idx(0);
		std::size_t r_idx(0);
		std::vector<std::size_t> on_split_points;
		on_split_points.reserve( 2 );
		for (const auto& idx : indices_) {
			auto p = points_[idx];
			if (l_idx < left_indices->size() && p[split_dim] < split) {
				(*left_indices)[l_idx++] = idx;
			}
			else if ( p[split_dim] > split ){
				(*right_indices)[r_idx++] = idx;
			}
			else{
				on_split_points.push_back( idx );
			}
		}
		assert( left_indices->size() + right_indices->size() == l_idx + r_idx + on_split_points.size() );

		for ( const auto& idx : on_split_points ) {
			if ( l_idx < left_indices->size() ) {
				(*left_indices)[l_idx++] = idx;
			}
			else {
				(*right_indices)[r_idx++] = idx;
			}
		}
		
		// todo: RAII?
		LeftChild* left_ptr = new decltype(left)(points_, *left_indices); //Used to circumvent the stack, which overflows
		left = *left_ptr;
		delete left_ptr;
		RightChild* right_ptr = new decltype(right)(points_, *right_indices);
		right = *right_ptr;
		delete right_ptr;
	}
	
	double split;

//	constexpr std::size_t n_left_points() { return n_points / 2; }
//	constexpr std::size_t n_right_points() { return n_points - (n_points / 2) };
	
	using LeftChild =
		typename std::conditional<(n_points / 2 <= max_points_per_node),
			KDTreeLeaf<CoordsType, dimension, n_points_total, n_points / 2, Me>,
			KDTree<CoordsType, dimension, n_points_total, n_points / 2, max_points_per_node, (split_dim + 1) % dimension, Me>
		>::type;
	using RightChild =
		typename std::conditional< (n_points - (n_points / 2) <= max_points_per_node),
			KDTreeLeaf<CoordsType, dimension, n_points_total, n_points - (n_points / 2), Me>,
			KDTree<CoordsType, dimension, n_points_total, n_points - (n_points / 2), max_points_per_node, (split_dim + 1) % dimension, Me>
		>::type;

	friend Parent;

	//<CoordsType, dimension, n_points / 2, max_points_per_node, (split_dim + 1) % dimension>
	//<CoordsType, dimension, n_points / 2, max_points_per_node, (split_dim + 1) % dimension>
	LeftChild left;
	RightChild right;

	void radius_search_ ( const Point& p, double radius, std::vector<std::size_t>& neighbors ) const  //TODO: make private
	{
		if ( p[split_dim] + radius < split ) {
			return left.radius_search_( p, radius, neighbors );
		}
		else if ( p[split_dim] - radius >= split ) {
			return right.radius_search_( p, radius, neighbors );
		}
		else {
			left.radius_search_( p, radius, neighbors );
			right.radius_search_( p, radius, neighbors );
			return;
		}
	}

};



}//namespace kdt

// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/Geometry
// Distributed under the MIT Software License (X11 license)
// (See accompanying file LICENSE)


#pragma once

#include <fplus/fplus.hpp>

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>

struct bgr_col
{
    bgr_col(std::uint8_t b, std::uint8_t g, std::uint8_t r) :
        b_(b), g_(g), r_(r)
    {}
    std::uint8_t b_;
    std::uint8_t g_;
    std::uint8_t r_;
};

inline bool operator == ( const bgr_col& lhs, const bgr_col& rhs ) {
	return (lhs.b_ == rhs.b_ && lhs.g_ == rhs.g_ && lhs.r_ == rhs.r_);
}

inline bool operator != ( const bgr_col& lhs, const bgr_col& rhs ) {
	return !(lhs == rhs);
}

inline bgr_col operator + ( const bgr_col& lhs, const bgr_col& rhs ) {
	return bgr_col(lhs.b_+rhs.b_, lhs.g_+rhs.g_, lhs.r_+rhs.r_);
}

inline bgr_col operator * ( const bgr_col&col, const double factor ) {
	std::uint8_t b = fplus::round<double,std::uint8_t>(static_cast<double>(col.b_) * factor);
	std::uint8_t g = fplus::round<double, std::uint8_t>( static_cast<double>(col.g_) * factor );
	std::uint8_t r = fplus::round<double, std::uint8_t>( static_cast<double>(col.r_) * factor );
	return { b,g,r };
}


struct img_pos
{
    img_pos(std::size_t x, std::size_t y) :
        x_(x), y_(y)
    {}
    std::size_t x_;
    std::size_t y_;
};

struct size_2d
{
public:
    size_2d(std::size_t width, std::size_t height) :
        width_(width), height_(height)
    {}
    std::size_t area() const { return width_ * height_; }
    std::size_t width_;
    std::size_t height_;
};

class bgr_image
{
public:
    bgr_image(const size_2d& size) :
        size_(size),
        data_(size_.area(), bgr_col(0, 0, 0))
    {};
	bgr_image( const size_2d& size, const bgr_col& col ) :
		size_( size ),
		data_( size_.area(), col )
	{};
    const size_2d& size() const { return size_; }
    const bgr_col& pix(const img_pos& pos) const {
		if ( pos.y_ >= size().height_ || pos.x_ >= size().width_ ) {
			throw std::out_of_range( "img_pos " + std::to_string( pos.x_ ) + "," + std::to_string( pos.x_ ) + " out of range!" );
		}
        return data_[pos.y_ * size().width_ + pos.x_];
    }
	const std::vector<bgr_col>& get_data() const{
		return data_;
	}
    bgr_col& pix(const img_pos& pos) {
		if ( pos.y_ >= size().height_ || pos.x_ >= size().width_ ) {
			throw std::out_of_range( "img_pos " + std::to_string(pos.x_) + "," + std::to_string(pos.y_) + " out of range!" );
		}
        return data_[pos.y_ * size().width_ + pos.x_];
    }
	void append_rows( const bgr_image& img ) {
		assert( img.size().width_ == size().width_ );
		size_.height_ += img.size().height_;
		const auto& add_data = img.get_data();
		data_ = fplus::append( data_, add_data );
	}
    bool save(const std::string& filepath) const;
private:
    size_2d size_;
    std::vector<bgr_col> data_;
};

inline bool bgr_image::save(const std::string& filepath) const
{
    std::ofstream file(filepath +".ppm", std::fstream::binary);
    if (file.bad())
        return false;
    file
        << "P6 "
        << std::to_string(size().width_)
        << " "
        << std::to_string(size().height_)
        << " "
        << "255"
<< "\n";
for ( std::size_t y = 0; y < size().height_; ++y )
{
	for ( std::size_t x = 0; x < size().width_; ++x )
	{
		const auto col = pix( img_pos( x, y ) );
		file << col.r_;
		file << col.g_;
		file << col.b_;
	}
}
return true;
}

namespace internal {
	inline std::vector<img_pos> line_pixel( const img_pos& p1, const img_pos& p2 ) {
		double dist = std::ceil( std::sqrt( (p2.x_ - p1.x_)*(p2.x_ - p1.x_) + (p2.y_ - p1.y_)*(p2.y_ - p1.y_) ) );

		std::vector<img_pos> result;
		result.reserve( static_cast<int>(dist) + 2 );

		int x1 = p1.x_; int y1 = p1.y_;
		int x2 = p2.x_; int y2 = p2.y_;

		int delta_x( x2 - x1 );
		// if x1 == x2, then it does not matter what we set here
		signed char const ix( (delta_x > 0) - (delta_x < 0) );
		delta_x = std::abs( delta_x ) << 1;

		int delta_y( y2 - y1 );
		// if y1 == y2, then it does not matter what we set here
		signed char const iy( (delta_y > 0) - (delta_y < 0) );
		delta_y = std::abs( delta_y ) << 1;

		result.push_back( img_pos( x1, y1 ) );

		if ( delta_x >= delta_y )
		{
			// error may go below zero
			int error( delta_y - (delta_x >> 1) );

			while ( x1 != x2 )
			{
				if ( (error >= 0) && (error || (ix > 0)) )
				{
					error -= delta_x;
					y1 += iy;
				}
				// else do nothing

				error += delta_y;
				x1 += ix;

				result.push_back( img_pos( x1, y1 ) );
			}
		}
		else
		{
			// error may go below zero
			int error( delta_x - (delta_y >> 1) );

			while ( y1 != y2 )
			{
				if ( (error >= 0) && (error || (iy > 0)) )
				{
					error -= delta_y;
					x1 += ix;
				}
				// else do nothing

				error += delta_x;
				y1 += iy;

				result.push_back( img_pos( x1, y1 ) );
			}
		}

		return result;
	}

	inline double line_dist( const img_pos& p1, const img_pos& p2, const img_pos& p ) {
		double x1 = static_cast<double>(p1.x_);
		double y1 = static_cast<double>(p1.y_);
		double x2 = static_cast<double>(p2.x_);
		double y2 = static_cast<double>(p2.y_);
		double x0 = static_cast<double>(p.x_);
		double y0 = static_cast<double>(p.y_);
		double nom = std::abs( (x2 - x1)*(y1 - y0) - (x1 - x0)*(y2 - y1) );
		double denom = std::sqrt( (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) );
		return nom / denom;
	}
}//namespace internal


inline void plot_pixel( bgr_image& image, const img_pos& p, const bgr_col& col ) {
	if ( p.y_ >= image.size().height_ || p.x_ >= image.size().width_ ) {
		return;
	}
	image.pix( p ) = col;
}

inline void plot_line_segment( bgr_image& image, const img_pos& p1, const img_pos& p2, const bgr_col& col )
{
	auto pixel = internal::line_pixel( p1, p2 );
	for ( const auto& p : pixel ) {
		plot_pixel( image, p, col );
	}
}

inline void plot_line_segment_antialiased( bgr_image& image, const img_pos& p1, const img_pos& p2, const bgr_col& col )
{
	auto pixel = internal::line_pixel( p1, p2 );
	for ( std::size_t idx = 1; idx < pixel.size()-1; idx+=2){
		auto pos = pixel[idx];
		for ( int y = -1; y <= 1; y++ ) {
			for ( int x = -1; x <= 1; x++ ) {
				auto p = img_pos( static_cast<std::size_t>(static_cast<int>(pos.x_) + x),
								  static_cast<std::size_t>(static_cast<int>(pos.y_) + y) );
				double dist = internal::line_dist( p1, p2, p );
				assert( dist >= 0.0 );
				if ( dist >= 0.99999 ) { continue; }
				bgr_col draw_col = image.pix(p)*dist + col*(1.0-dist);
				plot_pixel( image, p, draw_col );
			}
		}
	}
}

inline void plot_circle( bgr_image& image, const img_pos& center, std::size_t radius, const bgr_col& col)
{
	if ( radius == 0 ) {
		plot_pixel( image, center, col );
	}
    int x = static_cast<int>(radius);
    int y = 0;
    int x0 = center.x_;
    int y0 = center.y_;
    int err = 0;

    while (x >= y)
    {
        plot_pixel( image, img_pos(x0 + x, y0 + y), col );
        plot_pixel( image, img_pos(x0 + y, y0 + x), col );
        plot_pixel( image, img_pos(x0 - y, y0 + x), col );
        plot_pixel( image, img_pos(x0 - x, y0 + y), col );
        plot_pixel( image, img_pos(x0 - x, y0 - y), col );
        plot_pixel( image, img_pos(x0 - y, y0 - x), col );
        plot_pixel( image, img_pos(x0 + y, y0 - x), col );
        plot_pixel( image, img_pos(x0 + x, y0 - y), col );

        y += 1;
        if (err <= 0)
        {
            err += 2*y + 1;
        }
        if (err > 0)
        {
            x -= 1;
            err -= 2*x + 1;
        }
    }

}

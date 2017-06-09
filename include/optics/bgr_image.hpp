#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

struct bgr_col
{
    bgr_col(std::uint8_t b, std::uint8_t g, std::uint8_t r) :
        b_(b), g_(g), r_(r)
    {}
    std::uint8_t b_;
    std::uint8_t g_;
    std::uint8_t r_;
};

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
        data_(size_.area() * 3, bgr_col(0, 0, 0))
    {
    };
	bgr_image( const size_2d& size, const bgr_col& col ) :
		size_( size ),
		data_( size_.area() * 3, col )
	{
	};
    const size_2d& size() const { return size_; }
    const bgr_col& pix(const img_pos& pos) const {
		if ( pos.y_ >= size().height_ || pos.x_ >= size().width_ ) {
			throw std::out_of_range( "img_pos " + std::to_string( pos.x_ ) + "," + std::to_string( pos.x_ ) + " out of range!" );
		}
        return data_[pos.y_ * size().width_ + pos.x_];
    }
    bgr_col& pix(const img_pos& pos) {
		if ( pos.y_ >= size().height_ || pos.x_ >= size().width_ ) {
			throw std::out_of_range( "img_pos " + std::to_string(pos.x_) + "," + std::to_string(pos.y_) + " out of range!" );
		}
        return data_[pos.y_ * size().width_ + pos.x_];
    }
    bool save(const std::string& filepath) const;
private:
    size_2d size_;
    std::vector<bgr_col> data_;
};

bool bgr_image::save(const std::string& filepath) const
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
    for (std::size_t y = 0; y < size().height_; ++y)
    {
        for (std::size_t x = 0; x < size().width_; ++x)
        {
            const auto col = pix(img_pos(x, y));
            file << col.r_;
            file << col.g_;
            file << col.b_;
        }
    }
    return true;
}

void draw_pixel( bgr_image& image, const img_pos& p, const bgr_col& col ) {
	if ( p.y_ >= image.size().height_ || p.x_ >= image.size().width_ ) {
		return;
	}
	image.pix( p ) = col;
}

void plot_line_segment( bgr_image& image, const img_pos& p1, const img_pos& p2, const bgr_col& col )
{
	std::vector<img_pos> result;
	double dist = std::ceil( std::sqrt( (p2.x_ - p1.x_)*(p2.x_ - p1.x_) + (p2.y_ - p1.y_)*(p2.y_ - p1.y_)) );
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

	for ( const auto& p : result ) {
		draw_pixel( image, p, col );
	}
}

#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

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
    const size_2d& size() const { return size_; }
    const bgr_col& pix(const img_pos& pos) const {
        return data_[pos.y_ * size().width_ + pos.x_];
    }
    bgr_col& pix(const img_pos& pos) {
        const auto idx = pos.y_ * size().width_ + pos.x_;
        return data_[pos.y_ * size().width_ + pos.x_];
    }
    bool save(const std::string& filepath) const;
private:
    size_2d size_;
    std::vector<bgr_col> data_;
};

bool bgr_image::save(const std::string& filepath) const
{
    std::ofstream file(filepath, std::fstream::binary);
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

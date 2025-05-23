/*
 * easy_image.cc
 * Copyright (C) 2011  Daniel van den Akker
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "easy_image.h"
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <sstream>

#ifndef le32toh
#define le32toh(x) (x)
#endif

namespace
{
	//structs borrowed from wikipedia's article on the BMP file format
	struct bmpfile_magic
	{
			uint8_t magic[2];
	};

	struct bmpfile_header
	{
			uint32_t file_size;
			uint16_t reserved_1;
			uint16_t reserved_2;
			uint32_t bmp_offset;
	};
	struct bmp_header
	{
			uint32_t header_size;
			int32_t width;
			int32_t height;
			uint16_t nplanes;
			uint16_t bits_per_pixel;
			uint32_t compress_type;
			uint32_t pixel_size;
			int32_t hres;
			int32_t vres;
			uint32_t ncolors;
			uint32_t nimpcolors;
	};
	//copy-pasted from lparser.cc to allow these classes to be used independently from each other
	class enable_exceptions
	{
		private:
			std::ios& ios;
			std::ios::iostate state;
		public:
			enable_exceptions(std::ios& an_ios, std::ios::iostate exceptions) :
				ios(an_ios)
			{
				state = ios.exceptions();
				ios.exceptions(exceptions);
			}
			~enable_exceptions()
			{
				ios.exceptions(state);
			}
	};
	//helper function to convert a number (char, int, ...) to little endian
	//regardless of the endiannes of the system
	//more efficient machine-dependent functions exist, but this one is more portable
	template<typename T> T to_little_endian(T value)
	{
		//yes, unions must be used with caution, but this is a case in which a union is needed
		union
		{
				T t;
				uint8_t bytes[sizeof(T)];
		} temp_storage;

		for (uint8_t i = 0; i < sizeof(T); i++)
		{
			temp_storage.bytes[i] = value & 0xFF;
			value >>= 8;
		}
		return temp_storage.t;
	}

	template<typename T> T from_little_endian(T value)
	{
		//yes, unions must be used with caution, but this is a case in which a union is needed
		union
		{
				T t;
				uint8_t bytes[sizeof(T)];
		} temp_storage;
		temp_storage.t = value;
		T retVal = 0;

		for (uint8_t i = 0; i < sizeof(T); i++)
		{
			retVal = (retVal << 8) | temp_storage.bytes[sizeof(T) - i - 1];
		}
		return retVal;
	}

}
img::Color::Color() :
	blue(0), green(0), red(0)
{
}
img::Color::Color(uint8_t r, uint8_t g, uint8_t b) :
	blue(b), green(g), red(r)
{
}
img::Color::~Color()
{
}

img::UnsupportedFileTypeException::UnsupportedFileTypeException(std::string const& msg) :
	message(msg)
{
}
img::UnsupportedFileTypeException::UnsupportedFileTypeException(const UnsupportedFileTypeException &original)
: std::exception(original)
, message(original.message)
{
}
img::UnsupportedFileTypeException::~UnsupportedFileTypeException() throw ()
{
}
img::UnsupportedFileTypeException& img::UnsupportedFileTypeException::operator=(UnsupportedFileTypeException const& original)
{
	this->message = original.message;
	return *this;
}
const char* img::UnsupportedFileTypeException::what() const throw ()
{
	return message.c_str();
}

img::EasyImage::EasyImage() :
	width(0), height(0), bitmap()
{
}

img::EasyImage::EasyImage(unsigned int _width, unsigned int _height, Color color) :
	width(_width), height(_height), bitmap(width * height, color)
{
}

img::EasyImage::EasyImage(EasyImage const& img) :
	width(img.width), height(img.height), bitmap(img.bitmap)
{
}

img::EasyImage::~EasyImage()
{
	bitmap.clear();
}

img::EasyImage& img::EasyImage::operator=(img::EasyImage const& img)
{
	width = img.width;
	height = img.height;
	bitmap.assign(img.bitmap.begin(),img.bitmap.end());
	return (*this);
}

unsigned int img::EasyImage::get_width() const
{
	return width;
}

unsigned int img::EasyImage::get_height() const
{
	return height;
}

void img::EasyImage::clear(Color color)
{
	for (std::vector<Color>::iterator i = bitmap.begin(); i != bitmap.end(); i++)
	{
		*i = color;
	}
}

img::Color& img::EasyImage::operator()(unsigned int x, unsigned int y)
{
	assert(x < this->width);
	assert(y < this->height);
	return bitmap.at(x * height + y);
}

img::Color const& img::EasyImage::operator()(unsigned int x, unsigned int y) const
{
	assert(x < this->width);
	assert(y < this->height);
	return bitmap.at(x * height + y);
}
void img::EasyImage::draw_line(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, Color color)
{
	if (x0 >= this->width || y0 >= this->height || x1 >= this->width || y1 > this->height) {
		std::stringstream ss;
		ss << "Drawing line from (" << x0 << "," << y0 << ") to (" << x1 << "," << y1 << ") in image of width "
			<< this->width << " and height " << this->height;
		throw std::runtime_error(ss.str());
	}
	if (x0 == x1)
	{
		// Special case for vertical line
		for (int i = (int) (std::max(y0, y1) - std::min(y0, y1)); i >= 0; i--)
		{
			(*this)(x0, std::min(y0, y1) + i) = color;
		}
	}
	else if (y0 == y1)
	{
		// Special case for horizontal line
		for (int i = (int) (std::max(x0, x1) - std::min(x0, x1)); i >= 0; i--)
		{
			(*this)(std::min(x0, x1) + i, y0) = color;
		}
	}
	else
	{
		if (x0 > x1)
		{
			// Flip points if x1 > x0 to ensure x0 has the lowest value
			std::swap(x0, x1);
			std::swap(y0, y1);
		}
		double m = ((double)y1 - (double)y0) / ((double)x1 - (double)x0);
		if (-1.0 <= m && m <= 1.0)
		{
			for (int i = (int) (x1 - x0); i >= 0; i--)
			{
				(*this)(x0 + i, (unsigned int)round(y0 + m * i)) = color;
			}
		}
		else if (m > 1.0)
		{
			for (int i = (int) (y1 - y0); i >= 0; i--)
			{
				(*this)((unsigned int)round(x0 + (i / m)), y0 + i) = color;
			}
		}
		else if (m < -1.0)
		{
			for (int i = (int) (y0 - y1); i >= 0; i--)
			{
				(*this)((unsigned int)round(x0 - (i / m)), y0 - i) = color;
			}
		}
	}
}
void img::EasyImage::draw_zbuf_line(ZBuffer &buffer, EasyImage &image,unsigned int x0, unsigned int y0, double z0, unsigned int x1, unsigned int y1, double z1, const Color &color)
{
	if (x0 >= this->width || y0 >= this->height || x1 >= this->width || y1 > this->height) {
		std::stringstream ss;
		ss << "Drawing line from (" << x0 << "," << y0 << ") to (" << x1 << "," << y1 << ") in image of width "
			<< this->width << " and height " << this->height;
		throw std::runtime_error(ss.str());
	}

	if (x0 == x1 && y0 == y1) {
		double z = (1.0 / z0 + 1.0 / z1) / 2.0;
		if (z < buffer[x0][y0]) {
			image(x0, y0) = color;
			buffer[x0][y0] = z;
		}
	}
	else if (x0 == x1 || y0 == y1) {
		unsigned int dx, dy;
		if (x1 > x0) {
			dx = x1 - x0;
		} else {
			dx = x0 - x1;
		}

		if (y1 > y0) {
			dy = y1 - y0;
		} else {
			dy = y0 - y1;
		}

		unsigned int steps;

		if (dx == 0) {
			steps = dy;
		} else {
			steps = dx;
		}

		for (unsigned int i = 0; i <= steps; ++i) {
			unsigned int x;
			if (dx == 0) {
				x = x0;

			} else {
				if (x0 < x1) {
					x = x0 + i;
				} else {
					x = x0 - i;
				}
			}

			unsigned int y;
			if (dy == 0) {
				y = y0;
			} else {
				if (y0 < y1) {
					y = y0 + i;
				} else {
					y = y0 - i;
				}
			}

			double alpha = 1.0 - ((double)i / steps);
			double z = alpha / z0 + (1.0 - alpha) / z1;

			if (z < buffer[x][y]) {
				image(x, y) = color;
				buffer[x][y] = z;
			}
		}
	}
	else
	{
		if (x0 > x1)
		{
			// Flip points if x1 > x0 to ensure x0 has the lowest value
			std::swap(x0, x1);
			std::swap(y0, y1);
			std::swap(z0, z1);
		}
		double m = ((double)y1 - (double)y0) / ((double)x1 - (double)x0);
		if (-1.0 <= m && m <= 1.0)
		{
			double a = x1 - x0;
			for (unsigned int i = 0; i <= x1 - x0; i++)
			{
				double t = i / a;
				unsigned int x = x0 + i;
				unsigned int y = (unsigned int)round(y0 + m * i);
				double current = buffer[x][y];
				double z = (1 - t) / z0 + t / z1;

				if (z < current) {
					image(x, y) = color;
					buffer[x][y] = z;
				}
			}
		}
		else if (m > 1.0)
		{
			double a = y1 - y0;
			for (unsigned int i = 0; i <= y1 - y0; i++)
			{
				double t = i / a;
				unsigned int y = y0 + i;
				unsigned int x = (unsigned int)round(x0 + i / m);
				double current = buffer[x][y];
				double z = (1 - t) / z0 + t / z1;

				if (z < current) {
					image(x, y) = color;
					buffer[x][y] = z;
				}
			}
		}
		else if (m < -1.0)
		{
			double a = y0 - y1;
			for (unsigned int i = 0; i <= y0 - y1; i++)
			{
				double t = i / a;
				unsigned int y = y0 - i;
				unsigned int x = (unsigned int)round(x0 - i / m);
				double current = buffer[x][y];
				double z = (1 - t) / z0 + t / z1;

				if (z < current) {
					image(x, y) = color;
					buffer[x][y] = z;
				}
			}
		}
	}
}
void img::EasyImage::draw_zbuf_triag(ZBuffer &buffer, EasyImage &image, Vector3D const& A, Vector3D const& B, Vector3D const& C, double d, double dx, double dy, Color color) {
	Vector3D a;
	a.x = d*A.x/-A.z+dx;
	a.y = d*A.y/-A.z+dy;
	Vector3D b;
	b.x = d*B.x/-B.z+dx;
	b.y = d*B.y/-B.z+dy;
	Vector3D c;
	c.x = d*C.x/-C.z+dx;
	c.y = d*C.y/-C.z+dy;

	int yMin = lround(std::min(std::min(a.y, b.y),c.y)+0.5);
	int yMax = lround(std::max(std::max(a.y, b.y),c.y)-0.5);

	for (int yi = yMin; yi <= yMax; yi++) {
		double xlAB = std::numeric_limits<double>::infinity();
		double xlAC = std::numeric_limits<double>::infinity();
		double xlBC = std::numeric_limits<double>::infinity();
		double xrAB = -std::numeric_limits<double>::infinity();
		double xrAC = -std::numeric_limits<double>::infinity();
		double xrBC = -std::numeric_limits<double>::infinity();
		if ((yi-a.y)*(yi-b.y)<=0 && a.y != b.y) {
			double xi = b.x + (a.x-b.x)*((yi-b.y)/(a.y-b.y));
			xlAB = xrAB = xi;
		}
		if ((yi-a.y)*(yi-c.y)<=0 && a.y != c.y) {
			double xi = c.x + (a.x-c.x)*((yi-c.y)/(a.y-c.y));
			xlAC = xrAC = xi;
		}
		if ((yi-b.y)*(yi-c.y)<=0 && b.y != c.y) {
			double xi = c.x + (b.x-c.x)*((yi-c.y)/(b.y-c.y));
			xlBC = xrBC = xi;
		}
		int xl = lround(std::min(std::min(xlAB,xlAC),xlBC)+0.5);
		int xr = lround(std::max(std::max(xrAB,xrAC),xrBC)-0.5);

		double xG = (a.x+b.x+c.x)/3;
		double yG = (a.y+b.y+c.y)/3;

		double _1zG = 1/(3*A.z) + 1/(3*B.z) + 1/(3*C.z);

		Vector3D u = Vector3D::vector(B.x-A.x,B.y-A.y,B.z-A.z);
		Vector3D v = Vector3D::vector(C.x-A.x,C.y-A.y,C.z-A.z);
		double w1 = u.y*v.z - u.z*v.y;
		double w2 = u.z*v.x - u.x*v.z;
		double w3 = u.x*v.y - u.y*v.x;
		const double k = w1*A.x + w2*A.y + w3*A.z;
		const double dzdx = w1/(-d*k);
		const double dzdy = w2/(-d*k);

		for (int i = xl; i <= xr; i++) {
			double current = buffer[i][yi];
			double z = 1.0001*_1zG+(i-xG)*dzdx + (yi-yG)*dzdy;
			if (z < current) {
				image(i, yi) = color;
				buffer[i][yi] = z;
			}
		}
	}
}


std::ostream& img::operator<<(std::ostream& out, EasyImage const& image)
{

	//temporaryily enable exceptions on output stream
	enable_exceptions(out, std::ios::badbit | std::ios::failbit);
	//declare some struct-vars we're going to need:
	bmpfile_magic magic;
	bmpfile_header file_header;
	bmp_header header;
	uint8_t padding[] =
	{ 0, 0, 0, 0 };
	//calculate the total size of the pixel data
	unsigned int line_width = image.get_width() * 3; //3 bytes per pixel
	unsigned int line_padding = 0;
	if (line_width % 4 != 0)
	{
		line_padding = 4 - (line_width % 4);
	}
	//lines must be aligned to a multiple of 4 bytes
	line_width += line_padding;
	unsigned int pixel_size = image.get_height() * line_width;

	//start filling the headers
	magic.magic[0] = 'B';
	magic.magic[1] = 'M';

	file_header.file_size = to_little_endian(pixel_size + sizeof(file_header) + sizeof(header) + sizeof(magic));
	file_header.bmp_offset = to_little_endian(sizeof(file_header) + sizeof(header) + sizeof(magic));
	file_header.reserved_1 = 0;
	file_header.reserved_2 = 0;
	header.header_size = to_little_endian(sizeof(header));
	header.width = to_little_endian(image.get_width());
	header.height = to_little_endian(image.get_height());
	header.nplanes = to_little_endian(1);
	header.bits_per_pixel = to_little_endian(24);//3bytes or 24 bits per pixel
	header.compress_type = 0; //no compression
	header.pixel_size = pixel_size;
	header.hres = to_little_endian(11811); //11811 pixels/meter or 300dpi
	header.vres = to_little_endian(11811); //11811 pixels/meter or 300dpi
	header.ncolors = 0; //no color palette
	header.nimpcolors = 0;//no important colors

	//okay that should be all the header stuff: let's write it to the stream
	out.write((char*) &magic, sizeof(magic));
	out.write((char*) &file_header, sizeof(file_header));
	out.write((char*) &header, sizeof(header));

	//okay let's write the pixels themselves:
	//they are arranged left->right, bottom->top, b,g,r
	for (unsigned int i = 0; i < image.get_height(); i++)
	{
		//loop over all lines
		for (unsigned int j = 0; j < image.get_width(); j++)
		{
			//loop over all pixels in a line
			//we cast &color to char*. since the color fields are ordered blue,green,red they should be written automatically
			//in the right order
			out.write((char*) &image(j, i), 3 * sizeof(uint8_t));
		}
		if (line_padding > 0)
			out.write((char*) padding, line_padding);
	}
	//okay we should be done
	return out;
}
std::istream& img::operator>>(std::istream& in, EasyImage & image)
{
	enable_exceptions(in, std::ios::badbit | std::ios::failbit);
	//declare some struct-vars we're going to need
	bmpfile_magic magic;
	bmpfile_header file_header;
	bmp_header header;
	//a temp buffer for reading the padding at the end of each line
	uint8_t padding[] =
	{ 0, 0, 0, 0 };

	//read the headers && do some sanity checks
	in.read((char*) &magic, sizeof(magic));
	if (magic.magic[0] != 'B' || magic.magic[1] != 'M')
		throw UnsupportedFileTypeException("Could not parse BMP File: invalid magic header");
	in.read((char*) &file_header, sizeof(file_header));
	in.read((char*) &header, sizeof(header));
	if (le32toh(header.pixel_size) + le32toh(file_header.bmp_offset) != le32toh(file_header.file_size))
		throw UnsupportedFileTypeException("Could not parse BMP File: file size mismatch");
	if (le32toh(header.header_size) != sizeof(header))
		throw UnsupportedFileTypeException("Could not parse BMP File: Unsupported BITMAPV5HEADER size");
	if (le32toh(header.compress_type) != 0)
		throw UnsupportedFileTypeException("Could not parse BMP File: Only uncompressed BMP files can be parsed");
	if (le32toh(header.nplanes) != 1)
		throw UnsupportedFileTypeException("Could not parse BMP File: Only one plane should exist in the BMP file");
	if (le32toh(header.bits_per_pixel) != 24)
		throw UnsupportedFileTypeException("Could not parse BMP File: Only 24bit/pixel BMP's are supported");
	//if height<0 -> read top to bottom instead of bottom to top
	bool invertedLines = from_little_endian(header.height) < 0;
	image.height = std::abs(from_little_endian(header.height));
	image.width = std::abs(from_little_endian(header.width));
	unsigned int line_padding = from_little_endian(header.pixel_size) / image.height - (3 * image.width);
	//re-initialize the image bitmap
	image.bitmap.clear();
	image.bitmap.assign(image.height * image.width, Color());
	//okay let's read the pixels themselves:
	//they are arranged left->right., bottom->top if height>0, top->bottom if height<0, b,g,r
	for (unsigned int i = 0; i < image.get_height(); i++)
	{
		//loop over all lines
		for (unsigned int j = 0; j < image.get_width(); j++)
		{
			//loop over all pixels in a line
			//we cast &color to char*. since the color fields are ordered blue,green,red, the data read should be written in the right variables
			if (invertedLines)
			{
				//store top-to-bottom
				in.read((char*) &image(j, image.height - 1 - i), 3 * sizeof(uint8_t));
			}
			else
			{
				//store bottom-to-top
				in.read((char*) &image(j, i), 3 * sizeof(uint8_t));
			}
		}
		if (line_padding > 0)
		{
			in.read((char*) padding, line_padding);
		}
	}
	//okay we're done
	return in;
}

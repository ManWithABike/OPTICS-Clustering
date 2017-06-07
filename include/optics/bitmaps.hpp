// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)


#pragma once


#include <Windows.h>
#include <algorithm>
#include <memory>

// Save the bitmap to a bmp file  
void SaveBitmapToFile( BYTE* pBitmapBits,
					   LONG lWidth,
					   LONG lHeight,
					   WORD wBitsPerPixel,
					   const unsigned long& padding_size,
					   LPCTSTR lpszFileName )
{
	// Some basic bitmap parameters  
	unsigned long headers_size = sizeof( BITMAPFILEHEADER ) +
		sizeof( BITMAPINFOHEADER );

	unsigned long pixel_data_size = lHeight * ((lWidth * (wBitsPerPixel / 8)) + padding_size);

	BITMAPINFOHEADER bmpInfoHeader = { 0 };

	// Set the size  
	bmpInfoHeader.biSize = sizeof( BITMAPINFOHEADER );

	// Bit count  
	bmpInfoHeader.biBitCount = wBitsPerPixel;

	// Use all colors  
	bmpInfoHeader.biClrImportant = 0;

	// Use as many colors according to bits per pixel  
	bmpInfoHeader.biClrUsed = 0;

	// Store as un Compressed  
	bmpInfoHeader.biCompression = BI_RGB;

	// Set the height in pixels  
	bmpInfoHeader.biHeight = lHeight;

	// Width of the Image in pixels  
	bmpInfoHeader.biWidth = lWidth;

	// Default number of planes  
	bmpInfoHeader.biPlanes = 1;

	// Calculate the image size in bytes  
	bmpInfoHeader.biSizeImage = pixel_data_size;

	BITMAPFILEHEADER bfh = { 0 };

	// This value should be values of BM letters i.e 0x4D42  
	// 0x4D = M 0×42 = B storing in reverse order to match with endian  
	bfh.bfType = 0x4D42;
	//bfh.bfType = 'B'+('M' << 8); 

	// <<8 used to shift ‘M’ to end  */  

	// Offset to the RGBQUAD  
	bfh.bfOffBits = headers_size;

	// Total size of image including size of headers  
	bfh.bfSize = headers_size + pixel_data_size;

	// Create the file in disk to write  
	HANDLE hFile = CreateFile( lpszFileName,
							   GENERIC_WRITE,
							   0,
							   NULL,
							   CREATE_ALWAYS,
							   FILE_ATTRIBUTE_NORMAL,
							   NULL );

	// Return if error opening file  
	if ( !hFile ) return;

	DWORD dwWritten = 0;

	// Write the File header  
	WriteFile( hFile,
			   &bfh,
			   sizeof( bfh ),
			   &dwWritten,
			   NULL );

	// Write the bitmap info header  
	WriteFile( hFile,
			   &bmpInfoHeader,
			   sizeof( bmpInfoHeader ),
			   &dwWritten,
			   NULL );

	// Write the RGB Data  
	WriteFile( hFile,
			   pBitmapBits,
			   bmpInfoHeader.biSizeImage,
			   &dwWritten,
			   NULL );

	// Close the file handle  
	CloseHandle( hFile );
}
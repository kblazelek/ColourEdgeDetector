#include "stdafx.h"
#include <algorithm>
#include <fstream>
#include <boost/gil/image.hpp>
#include <boost/gil/typedefs.hpp>
#include <boost/gil/color_convert.hpp>
#include <boost/gil/extension/io/jpeg_io.hpp>
#include <boost/math/quaternion.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::gil;
using namespace boost::math;
using namespace std;

vector<vector<quaternion<double>>> ConvertRGBImageToQuaternionMatrix(rgb8_image_t &image)
{
	// Obraz 768x512x3 jest zamieniany na macierz kwaternionow o rozmiarze 768x512.
	// Dla pixela o barwie (R,G,B) = (10,20,30) kwaternion (czysty) dla 24 bitowego obrazu wyglada nastepujaco:
	// Q = 0 + (10 - 127,5)b + (20 - 127,5)c + (30 - 127,5)d
	cout << "Zamiana obrazu RGB na macierz kwaternionow... " << "\n";
	int imageWidth = image.width();
	int imageHeight = image.height();
	quaternion<double> emptyQuaternion(0, 0, 0, 0);
	vector<quaternion<double>> rowOfEmptyQuaternions(imageWidth, emptyQuaternion);
	vector<vector<quaternion<double>>> quaternionMatrix(imageHeight, rowOfEmptyQuaternions);
	for (int row = 0; row < imageHeight; ++row)
	{
		rgb8_view_t::x_iterator rowIterator = image._view.row_begin(row);
		for (int column = 0; column < imageWidth; ++column)
		{
			quaternionMatrix[row][column] = quaternion<double>(0, rowIterator[column][0] - 127.5, rowIterator[column][1] - 127.5, rowIterator[column][2] - 127.5);
		}
	}
	return quaternionMatrix;
}

rgb8_image_t ConvertQuaternionMatrixToRGBImage(vector<vector<quaternion<double>>> quaternionMatrix)
{
	cout << "Zamiana macierzy kwaternionow na obraz RGB... " << "\n";
	int imageHeight = quaternionMatrix.size();
	int imageWidth = quaternionMatrix[0].size();
	rgb8_image_t img(imageWidth, imageHeight);
	rgb8_image_t::view_t v = view(img);
	for (int x = 0; x < imageWidth; ++x)
	{
		for (int y = 0; y < imageHeight; ++y)
		{
			double r = quaternionMatrix[y][x].R_component_2() + 127.5;
			double g = quaternionMatrix[y][x].R_component_3() + 127.5;
			double b = quaternionMatrix[y][x].R_component_4() + 127.5;
			v(x, y) = rgb8_pixel_t((unsigned char)r, (unsigned char)g, (unsigned char)b);
		}
	}
	return img;
}

quaternion<double> normalizeQuaternion(const quaternion<double> &q)
{
	double w = q.R_component_1();
	double x = q.R_component_2();
	double y = q.R_component_3();
	double z = q.R_component_4();
	double inv_len = 1.0f / std::sqrt(w * w + x * x + y * y + z * z);
	return q * inv_len;
}

vector<vector<quaternion<double>>> ApplyHypercomplexFilter(vector<vector<quaternion<double>>> quaternionMatrix, rgb8_pixel_t color1, rgb8_pixel_t color2)
{
	cout << "Stosowanie filtru hiperzespolonego... " << "\n";
	int imageHeight = quaternionMatrix.size();
	int imageWidth = quaternionMatrix[0].size();
	quaternion<double> emptyQuaternion(0, 0, 0, 0);
	vector<quaternion<double>> rowOfEmptyQuaternionsForMasks(3, emptyQuaternion);
	vector<quaternion<double>> rowOfEmptyQuaternionsForImage(imageWidth, emptyQuaternion);
	vector<vector<quaternion<double>>> leftMask(3, rowOfEmptyQuaternionsForMasks);
	vector<vector<quaternion<double>>> rightMask(3, rowOfEmptyQuaternionsForMasks);
	vector<vector<quaternion<double>>> imageAfterFilter(imageHeight, rowOfEmptyQuaternionsForImage);
	quaternion<double> quaternionC1(0, color1[0] - 127.5, color1[1] - 127.5, color1[2] - 127.5);
	quaternion<double> quaternionC2(0, color2[0] - 127.5, color2[1] - 127.5, color2[2] - 127.5);
	quaternion<double> normalizedQuaternionC1 = normalizeQuaternion(quaternionC1);
	quaternion<double> normalizedQuaternionC2 = normalizeQuaternion(quaternionC2);
	leftMask[0][0] = normalizedQuaternionC2 / sqrt(6);
	leftMask[0][1] = normalizedQuaternionC2 / sqrt(6);
	leftMask[0][2] = normalizedQuaternionC2 / sqrt(6);
	leftMask[2][0] = 1 / sqrt(6);
	leftMask[2][1] = 1 / sqrt(6);
	leftMask[2][2] = 1 / sqrt(6);
	rightMask[0][0] = 1 / sqrt(6);
	rightMask[0][1] = 1 / sqrt(6);
	rightMask[0][2] = 1 / sqrt(6);
	rightMask[2][0] = normalizedQuaternionC1 / sqrt(6);
	rightMask[2][1] = normalizedQuaternionC1 / sqrt(6);
	rightMask[2][2] = normalizedQuaternionC1 / sqrt(6);

	// Zaczynamy od 2 indeksu (3 pixela) od lewej i gory, bo nie policzymy maski dla wczeœniejszych (nie maja sasiadow na lewo/od gory)
	for (int s = 2; s < imageWidth; ++s)
	{
		for (int t = 2; t < imageHeight; ++t)
		{
			//cout << "Pixel " << s << "\t " << t << "\n";
			for (int x = 0; x < 3; ++x)
			{
				for (int y = 0; y < 3; ++y)
				{
					imageAfterFilter[t][s] += leftMask[x][y] * quaternionMatrix[t - y][s - x] * rightMask[x][y];
				}
			}
		}
	}
	return imageAfterFilter;
}
int main()
{
	rgb8_image_t img;
	jpeg_read_image("..\\..\\Images\\tulips.jpg", img);
	//jpeg_read_image("..\\..\\Images\\colours.jpg", img);
	int imageWidth = img.width();
	int imageHeight = img.height();
	printf("Image width: %i\n", imageWidth);
	printf("Image height: %i\n", imageHeight);
	vector<vector<quaternion<double>>> quaternionMatrix = ConvertRGBImageToQuaternionMatrix(img);
	rgb8_pixel_t color1(237, 157, 157); // czerwony
	rgb8_pixel_t color2(255, 255, 255); // bialy
	vector<vector<quaternion<double>>> quaternionMatrixAfterFilter = ApplyHypercomplexFilter(quaternionMatrix, color1, color2);
	rgb8_image_t imageAfterConversion = ConvertQuaternionMatrixToRGBImage(quaternionMatrixAfterFilter);
	jpeg_write_view("imageAfterHypercomplexFilter.jpg", view(imageAfterConversion));
	getchar();
	return 0;
}


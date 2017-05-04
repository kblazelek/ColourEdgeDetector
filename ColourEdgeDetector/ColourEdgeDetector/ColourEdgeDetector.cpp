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

gray8_image_t ConvertQuaternionMatrixToGrayImage(vector<vector<quaternion<double>>> quaternionMatrix)
{
	cout << "Zamiana macierzy kwaternionow na obraz w skali szarosci... " << "\n";
	int imageHeight = quaternionMatrix.size();
	int imageWidth = quaternionMatrix[0].size();
	gray8_image_t img(imageWidth, imageHeight);
	gray8_image_t::view_t v = view(img);
	for (int x = 0; x < imageWidth; ++x)
	{
		for (int y = 0; y < imageHeight; ++y)
		{
			double a = quaternionMatrix[y][x].R_component_1() + 127.5;
			v(x, y) = gray8_pixel_t((unsigned char)a);
		}
	}
	return img;
}

quaternion<double> normalizeQuaternion(const quaternion<double> &q)
{
	double a = q.R_component_1();
	double i = q.R_component_2();
	double j = q.R_component_3();
	double k = q.R_component_4();
	double inv_len = 1.0f / std::sqrt(a * a + i * i + j * j + k * k);
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
	leftMask[2][0] = 1.0f / sqrt(6);
	leftMask[2][1] = 1.0f / sqrt(6);
	leftMask[2][2] = 1.0f / sqrt(6);
	rightMask[0][0] = 1.0f / sqrt(6);
	rightMask[0][1] = 1.0f / sqrt(6);
	rightMask[0][2] = 1.0f / sqrt(6);
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
					imageAfterFilter[t][s] += leftMask[y][x] * quaternionMatrix[t - y][s - x] * rightMask[y][x];
				}
			}
		}
	}
	return imageAfterFilter;
}
// <KAMIL>
quaternion<double> convertColorToQuaternion(rgb8_pixel_t color) {
	quaternion<double> quaternion(0, 
								  color[0] - MID_GREY_COLOR, 
								  color[1] - MID_GREY_COLOR, 
								  color[2] - MID_GREY_COLOR);
	return quaternion;
}

double quaternionMagnitude(const quaternion<double> quaternion) {
	return sqrt(pow(quaternion.R_component_1(), 2)
			  + pow(quaternion.R_component_2(), 2)
			  + pow(quaternion.R_component_3(), 2)
			  + pow(quaternion.R_component_4(), 2));
}

double calculateQuaternionScalar(rgb8_pixel_t color1, rgb8_pixel_t color2) {
	return -(quaternionMagnitude(convertColorToQuaternion(color1))
		   + quaternionMagnitude(convertColorToQuaternion(color2)))
		   / 2;
}

bool compareQuaternionScalarCriteria(const quaternion<double> quaternion, rgb8_pixel_t color1, rgb8_pixel_t color2) {
	return calculateQuaternionScalar(color1, color2) == quaternion.R_component_1() ? true : false;
}

double quaternionVector(const quaternion<double> quaternion) {
	return sqrt(pow(quaternion.R_component_2(), 2)
			  + pow(quaternion.R_component_3(), 2)
			  + pow(quaternion.R_component_4(), 2));
}

double calculateQuaternionVector(const quaternion<double> quaternion) {
	return quaternionVector(quaternion);
}

bool compareQuaternionVectorCriteria(const quaternion<double> quaternion) {
	return calculateQuaternionVector(quaternion) == 0 ? true : false;
}

bool compareNormalizedQuaternionScalarCriteria(quaternion<double> quaternion) {
	return normalizeQuaternion(quaternion).R_component_1 < 0 ? true : false;
}

bool compareNormalizedNonRelaxedQuaternionVectorCriteria(const quaternion<double> quaternion) {
	return quaternionVector(normalizeQuaternion(quaternion)) == 0 ? true : false;
}

void detectIdealEdgeByThresholding() {
	if (compareQuaternionScalarCriteria && compareQuaternionVectorCriteria) {

	}
}
// </KAMIL>
int main()
{
	rgb8_image_t img;
	//jpeg_read_image("..\\..\\Images\\tulips.jpg", img);
	jpeg_read_image("..\\..\\Images\\colours.jpg", img);
	int imageWidth = img.width();
	int imageHeight = img.height();
	printf("Image width: %i\n", imageWidth);
	printf("Image height: %i\n", imageHeight);
	vector<vector<quaternion<double>>> quaternionMatrix = ConvertRGBImageToQuaternionMatrix(img);
	rgb8_pixel_t color1(255, 0, 0); // czerwony
	rgb8_pixel_t color2(0, 0, 255); // niebieski
	vector<vector<quaternion<double>>> quaternionMatrixAfterFilter = ApplyHypercomplexFilter(quaternionMatrix, color1, color2);
	rgb8_image_t imageAfterConversion = ConvertQuaternionMatrixToRGBImage(quaternionMatrixAfterFilter);
	jpeg_write_view("RGBImageAfterHypercomplexFilter.jpg", view(imageAfterConversion));
	gray8_image_t test = ConvertQuaternionMatrixToGrayImage(quaternionMatrixAfterFilter);
	jpeg_write_view("GrayScaleImageAfterHypercomplexFilter.jpg", view(test));
	cout << "Gotowe\n";
	getchar();
	return 0;
}


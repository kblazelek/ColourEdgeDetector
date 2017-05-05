#include "stdafx.h"
#include <algorithm>
#include <fstream>
#include <boost/gil/image.hpp>
#include <boost/gil/typedefs.hpp>
#include <boost/gil/color_convert.hpp>
#include <boost/gil/extension/io/jpeg_io.hpp>
#include <boost/math/quaternion.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>

using namespace boost::gil;
using namespace boost::math;
using namespace std;
const double MID_GREY_COLOR = 127.5;

quaternion<double> ConvertColorToQuaternion(rgb8_pixel_t color)
{
	quaternion<double> quaternion(0,
		color[0] - MID_GREY_COLOR,
		color[1] - MID_GREY_COLOR,
		color[2] - MID_GREY_COLOR);
	return quaternion;
}


rgb8_pixel_t ConvertQuaternionToColor(quaternion<double> quaternion)
{
	unsigned int red = (unsigned int)(quaternion.R_component_2() + MID_GREY_COLOR);
	unsigned int green = (unsigned int)(quaternion.R_component_3() + MID_GREY_COLOR);
	unsigned int blue = (unsigned int)(quaternion.R_component_4() + MID_GREY_COLOR);
	return rgb8_pixel_t(red, green, blue);
}

vector<vector<quaternion<double>>> ConvertRGBImageToQuaternionMatrix(rgb8_image_t &image)
{
	// Obraz 768x512x3 jest zamieniany na macierz kwaternionow o rozmiarze 768x512.
	// Dla pixela o barwie (R,G,B) = (10,20,30) kwaternion (czysty) dla 24 bitowego obrazu wyglada nastepujaco:
	// Q = 0 + (10 - 127,5)b + (20 - 127,5)c + (30 - 127,5)d
	// (127,5, 127,5, 127,5) to centrum skali szarosci
	cout << "Zamiana obrazu RGB na macierz kwaternionow... " << "\n";
	int imageWidth = image.width();
	int imageHeight = image.height();
	quaternion<double> emptyQuaternion(0, 0, 0, 0);

	// Tworzymy macierz pustych kwaternionow o wymiarze badanego obrazka
	vector<quaternion<double>> columnOfEmptyQuaternions(imageHeight, emptyQuaternion);
	vector<vector<quaternion<double>>> quaternionMatrix(imageWidth, columnOfEmptyQuaternions);

	// Kolumna po kolumnie zamieniamy RGB na czyste kwaterniony
	for (int column = 0; column < imageWidth; ++column)
	{
		rgb8_view_t::y_iterator columnIterator = image._view.col_begin(column);
		for (int row = 0; row < imageHeight; ++row)
		{
			quaternionMatrix[column][row] = ConvertColorToQuaternion(columnIterator[row]);
		}
	}
	return quaternionMatrix;
}

rgb8_image_t ConvertQuaternionMatrixToRGBImage(vector<vector<quaternion<double>>> quaternionMatrix)
{
	cout << "Zamiana macierzy kwaternionow na obraz RGB... " << "\n";
	int imageWidth = quaternionMatrix.size();
	int imageHeight = quaternionMatrix[0].size();
	rgb8_image_t img(imageWidth, imageHeight);
	rgb8_image_t::view_t v = view(img);
	for (int x = 0; x < imageWidth; ++x)
	{
		for (int y = 0; y < imageHeight; ++y)
		{
			v(x, y) = ConvertQuaternionToColor(quaternionMatrix[x][y]);
		}
	}
	return img;
}

gray8_image_t ConvertQuaternionMatrixToGrayImage(vector<vector<quaternion<double>>> quaternionMatrix)
{
	cout << "Zamiana macierzy kwaternionow na obraz w skali szarosci... " << "\n";
	int imageWidth = quaternionMatrix.size();
	int imageHeight = quaternionMatrix[0].size();
	gray8_image_t img(imageWidth, imageHeight);
	gray8_image_t::view_t v = view(img);
	for (int x = 0; x < imageWidth; ++x)
	{
		for (int y = 0; y < imageHeight; ++y)
		{
			double a = quaternionMatrix[x][y].R_component_1() + MID_GREY_COLOR;
			v(x, y) = gray8_pixel_t((unsigned char)a);
		}
	}
	return img;
}

quaternion<double> NormalizeQuaternion(const quaternion<double> &q)
{
	double a = q.R_component_1();
	double i = q.R_component_2();
	double j = q.R_component_3();
	double k = q.R_component_4();
	double inv_len = 1.0f / std::sqrt(a * a + i * i + j * j + k * k);
	return q * inv_len;
}
gray8_image_t ApplyRelaxedThreshold(vector<vector<quaternion<double>>> quaternionMatrix, rgb8_pixel_t color1, rgb8_pixel_t color2)
{
	cout << "Progowanie i zamiana macierzy kwaternionow na obraz w skali szarosci... " << "\n";
	int imageWidth = quaternionMatrix.size();
	int imageHeight = quaternionMatrix[0].size();
	gray8_image_t img(imageWidth, imageHeight);
	gray8_image_t::view_t v = view(img);
	// Zamiana na czyste kwaterniony kolorow, pomiedzy ktorymi szukamy krawedzi
	quaternion<double> quaternionC1 = ConvertColorToQuaternion(color1);
	quaternion<double> quaternionC2 = ConvertColorToQuaternion(color2);

	// Normalizacja czystych kwaternionow
	quaternion<double> normalizedQuaternionC1 = NormalizeQuaternion(quaternionC1);
	quaternion<double> normalizedQuaternionC2 = NormalizeQuaternion(quaternionC2);

	// Liczenie wartosci progowej
	double a1 = normalizedQuaternionC1.R_component_1();
	double i1 = normalizedQuaternionC1.R_component_2();
	double j1 = normalizedQuaternionC1.R_component_3();
	double k1 = normalizedQuaternionC1.R_component_4();
	double a2 = normalizedQuaternionC2.R_component_1();
	double i2 = normalizedQuaternionC2.R_component_2();
	double j2 = normalizedQuaternionC2.R_component_3();
	double k2 = normalizedQuaternionC2.R_component_4();
	double cosinusBetweenColours = a1*a2 + i1*i2 + j1*j2 + k1*k2;
	double angleBetweenColours = acos(cosinusBetweenColours);
	double thresholdValue = min(0.204, tan(angleBetweenColours / 2));

	// Progowanie
	for (int x = 0; x < imageWidth; ++x)
	{
		for (int y = 0; y < imageHeight; ++y)
		{
			double a = quaternionMatrix[x][y].R_component_1();
			double i = quaternionMatrix[x][y].R_component_2();
			double j = quaternionMatrix[x][y].R_component_3();
			double k = quaternionMatrix[x][y].R_component_4();
			double vectorModule = sqrt(i*i + j*j + k*k);
			double scalarModule = sqrt(a*a);
			if (a < 0 && vectorModule/scalarModule < thresholdValue) // Mamy krawedz
			{
				v(x, y) = gray8_pixel_t(0); // Czarny pixel
			}
			else // Nie mamy krawedzi
			{
				v(x, y) = gray8_pixel_t(255); // Bialy pixel
			}
		}
	}
	return img;
}
vector<vector<quaternion<double>>> ApplyHypercomplexFilter(vector<vector<quaternion<double>>> quaternionMatrix, rgb8_pixel_t color1, rgb8_pixel_t color2)
{
	cout << "Stosowanie filtru hiperzespolonego... " << "\n";
	int imageWidth = quaternionMatrix.size();
	int imageHeight = quaternionMatrix[0].size();
	quaternion<double> emptyQuaternion(0, 0, 0, 0);
	vector<quaternion<double>> columnOfEmptyQuaternionsForMasks(3, emptyQuaternion);
	vector<quaternion<double>> columnOfEmptyQuaternionsForImage(imageHeight, emptyQuaternion);
	vector<vector<quaternion<double>>> leftMask(3, columnOfEmptyQuaternionsForMasks);
	vector<vector<quaternion<double>>> rightMask(3, columnOfEmptyQuaternionsForMasks);
	vector<vector<quaternion<double>>> imageAfterFilter(imageWidth, columnOfEmptyQuaternionsForImage);

	// Zamiana na czyste kwaterniony kolorow, pomiedzy ktorymi szukamy krawedzi
	quaternion<double> quaternionC1 = ConvertColorToQuaternion(color1);
	quaternion<double> quaternionC2 = ConvertColorToQuaternion(color2);

	// Normalizacja czystych kwaternionow
	quaternion<double> normalizedQuaternionC1 = NormalizeQuaternion(quaternionC1);
	quaternion<double> normalizedQuaternionC2 = NormalizeQuaternion(quaternionC2);

	// Lewa maska
	leftMask[0][0] = normalizedQuaternionC2 / sqrt(6);
	leftMask[1][0] = normalizedQuaternionC2 / sqrt(6);
	leftMask[2][0] = normalizedQuaternionC2 / sqrt(6);
	leftMask[0][2] = 1.0f / sqrt(6);
	leftMask[1][2] = 1.0f / sqrt(6);
	leftMask[2][2] = 1.0f / sqrt(6);

	// Prawa maska
	rightMask[0][0] = 1.0f / sqrt(6);
	rightMask[1][0] = 1.0f / sqrt(6);
	rightMask[2][0] = 1.0f / sqrt(6);
	rightMask[0][2] = normalizedQuaternionC1 / sqrt(6);
	rightMask[1][2] = normalizedQuaternionC1 / sqrt(6);
	rightMask[2][2] = normalizedQuaternionC1 / sqrt(6);

	// Zaczynamy od 2 indeksu (3 pixela) od lewej i gory, bo nie policzymy maski dla wczesniejszych (nie maja sasiadow na lewo/od gory)
	for (int s = 2; s < imageWidth; ++s) // Dla kazdej kolumny macierzy kwaternionow
	{
		for (int t = 2; t < imageHeight; ++t) // Dla kazdego wiersza macierzy kwaternionow
		{
			//cout << "Pixel " << s << "\t " << t << "\n";
			for (int x = 0; x < 3; ++x) // Dla kazdej kolumny maski
			{
				for (int y = 0; y < 3; ++y) // Dla kazdego wiersza maski
				{
					imageAfterFilter[s][t] += leftMask[x][y] * quaternionMatrix[s - x][t - y] * rightMask[x][y];
				}
			}
		}
	}
	return imageAfterFilter;
}

int main()
{
	rgb8_image_t img;

	//jpeg_read_image("..\\..\\Images\\tulips.jpg", img);
	//rgb8_pixel_t color1(237, 157, 157); // odcien czerwonego
	//rgb8_pixel_t color2(255, 255, 255); // bialy

	jpeg_read_image("..\\..\\Images\\colours.jpg", img);
	rgb8_pixel_t color1(255, 0, 0); // czerwony
	rgb8_pixel_t color2(0, 0, 255); // niebieski

	int imageWidth = img.width();
	int imageHeight = img.height();
	printf("Image width: %i\n", imageWidth);
	printf("Image height: %i\n", imageHeight);
	vector<vector<quaternion<double>>> quaternionMatrix = ConvertRGBImageToQuaternionMatrix(img);

	vector<vector<quaternion<double>>> quaternionMatrixAfterFilter = ApplyHypercomplexFilter(quaternionMatrix, color1, color2);
	rgb8_image_t vectorImage = ConvertQuaternionMatrixToRGBImage(quaternionMatrixAfterFilter);
	jpeg_write_view("RGBImageAfterHypercomplexFilter.jpg", view(vectorImage));
	gray8_image_t scalarImage = ConvertQuaternionMatrixToGrayImage(quaternionMatrixAfterFilter);
	jpeg_write_view("GrayScaleImageAfterHypercomplexFilter.jpg", view(scalarImage));
	gray8_image_t edgeImage = ApplyRelaxedThreshold(quaternionMatrixAfterFilter, color1, color2);
	jpeg_write_view("Edge.jpg", view(edgeImage));
	cout << "Gotowe\n";
	getchar();
	return 0;
}


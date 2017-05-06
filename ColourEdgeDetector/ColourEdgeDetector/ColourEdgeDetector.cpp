#include "stdafx.h"
#include <algorithm>
#include <fstream>
#include <boost/gil/image.hpp>
#include <boost/gil/typedefs.hpp>
#include <boost/gil/color_convert.hpp>
#include <boost/gil/extension/io/jpeg_io.hpp>
#include <boost/math/quaternion.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "boost/filesystem.hpp"
#include <cmath>

using namespace boost::filesystem;
using namespace boost::gil;
using namespace boost::math;
using namespace std;
const double MID_GREY_COLOR = 127.5;
string inputImage;
string outputDir;
bool showScalarPart = false;
bool showVectorPart = false;
string colour1Red;
string colour1Green;
string colour1Blue;
string colour2Red;
string colour2Green;
string colour2Blue;
rgb8_pixel_t colour1;
rgb8_pixel_t colour2;

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

quaternion<double> NormalizeQuaternion(const quaternion<double> &q)
{
	double a = q.R_component_1();
	double i = q.R_component_2();
	double j = q.R_component_3();
	double k = q.R_component_4();
	double inv_len = 1.0f / std::sqrt(a * a + i * i + j * j + k * k);
	return q * inv_len;
}

vector<vector<quaternion<double>>> ConvertRGBImageToQuaternionMatrix(rgb8_image_t &image)
{
	// Obraz 768x512x3 jest zamieniany na macierz kwaternionow o rozmiarze 768x512.
	// Dla pixela o barwie (R,G,B) = (10,20,30) kwaternion (czysty) dla 24 bitowego obrazu wyglada nastepujaco:
	// Q = 0 + (10 - 127,5)b + (20 - 127,5)c + (30 - 127,5)d
	// (127,5, 127,5, 127,5) to centrum skali szarosci
	cout << "Converting RGB image to quaternion matrix... " << "\n";
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
	cout << "Converting vector part of quaternion matrix to RGB image... " << "\n";
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
	cout << "Converting scalar part of quaternion matrix to grayscale image... " << "\n";
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

gray8_image_t ApplyRelaxedThreshold(vector<vector<quaternion<double>>> quaternionMatrix, rgb8_pixel_t color1, rgb8_pixel_t color2)
{
	cout << "Thresholding and saving detected edges as grayscale image... " << "\n";
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
			if (a < 0 && vectorModule / scalarModule < thresholdValue) // Mamy krawedz
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
	cout << "Applying hypercomplex filter to detect horizontal edges... " << "\n";
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

void ShowHelp()
{
	std::cerr << "Options:\n"
		<< "\t--help\t\tShow this help message" << std::endl
		<< "\t--inputImage PathToInputImage\t\tMandatory. Path to input image" << std::endl
		<< "\t--colour1 red green blue\t\tMandatory. First color between which to look for edges" << std::endl
		<< "\t--colour2 red green blue\t\tMandatory. Second color between which to look for edges" << std::endl
		<< "\t--outputDir PathToOutputDir\t\tOptional. Path to output directory" << std::endl
		<< "\t--showScalarPart\t\tOptional. Save scalar part of quaternion matrix as grayscale image" << std::endl
		<< "\t--showVectorPart\t\tOptional. Save vector part of quaternion matrix as RGB image" << std::endl
		<< "Example 1:" << std::endl
		<< "--inputImage \"C:\\colours.jpg\" --colour1 255 0 0 --colour2 0 0 255 --showScalarPart --showVectorPart" << std::endl
		<< "Example 2:" << std::endl
		<< "--inputImage \"C:\\colours.jpg\" --colour1 255 0 0 --colour2 0 0 255 --showScalarPart --showVectorPart --outputDir \"C:\\Output\"" << std::endl
		<< "Remarks:" << std::endl
		<< "Only 24 bit jpg images are supported as input image." << std::endl;
}

bool GetParameters(int argc, char* argv[])
{
	if (argc < 11) {
		ShowHelp();
		return false;
	}
	for (int i = 1; i < argc; ++i)
	{
		std::string arg = argv[i];
		if (arg == "--help")
		{
			ShowHelp();
			return false;
		}
		else if (arg == "--inputImage")
		{
			if (i + 1 < argc)
			{
				inputImage = argv[++i];
			}
			else
			{
				std::cerr << "--inputImage option requires one argument" << std::endl;
				return false;
			}
		}
		else if (arg == "--outputDir")
		{
			if (i + 1 < argc)
			{
				outputDir = argv[++i];
			}
			else
			{
				std::cerr << "--outputDir option requires one argument" << std::endl;
				return false;
			}
		}
		else if (arg == "--colour1")
		{
			if (i + 3 < argc)
			{
				colour1Red = argv[++i];
				colour1Green = argv[++i];
				colour1Blue = argv[++i];
			}
			else
			{
				std::cerr << "--colour1 option requires three arguments" << std::endl;
				return false;
			}
		}
		else if (arg == "--colour2")
		{
			if (i + 3 < argc)
			{
				colour2Red = argv[++i];
				colour2Green = argv[++i];
				colour2Blue = argv[++i];
			}
			else
			{
				std::cerr << "--colour2 option requires three arguments" << std::endl;
				return false;
			}
		}
		else if (arg == "--showScalarPart")
		{
			showScalarPart = true;
		}
		else if (arg == "--showVectorPart")
		{
			showVectorPart = true;
		}
	}
	return true;
}

bool ValidateParameters()
{
	if (inputImage.empty())
	{
		std::cerr << "--inputImage option is mandatory" << std::endl;
		return false;
	}
	if (!exists(inputImage))
	{
		std::cerr << "Provided input image does not exist" << std::endl;
		return false;
	}
	if (extension(inputImage) != ".jpg")
	{
		std::cerr << "Provided input image is not a jpg image. Only 24 bit jpg images are supported" << std::endl;
		return false;
	}
	if (!outputDir.empty())
	{
		if (outputDir.back() != '\\')
		{
			outputDir.push_back('\\');
		}
		if (!is_directory(outputDir))
		{
			cout << "Directory " << outputDir << " does not exist. Creating it right now" << std::endl;
			try
			{
				create_directories(outputDir);
			}
			catch (exception& e)
			{
				cout << "Could not create directory " << outputDir << std::endl;
				return false;
			}
		}
	}
	unsigned long red;
	unsigned long green;
	unsigned long blue;

	// Validate colour 1
	try
	{
		red = stoul(colour1Red);
		green = stoul(colour1Green);
		blue = stoul(colour1Blue);
		if (red > 255 || green > 255 || blue > 255)
		{
			cout << "Wrong colour range for colour 1. Colour ranges must be from 0 to 255";
		}
		else
		{
			colour1 = rgb8_pixel_t(red, green, blue);
		}
	}
	catch(exception& e)
	{
		cout << "Could not parse colour 1. Please provide colour in valid format, for example --colour1 255 255 255";
		return false;
	}

	// Validate colour 2
	try
	{
		red = stoul(colour2Red);
		green = stoul(colour2Green);
		blue = stoul(colour2Blue);
		if (red > 255 || green > 255 || blue > 255)
		{
			cout << "Wrong colour range for colour 2. Colour ranges must be from 0 to 255";
		}
		else
		{
			colour2 = rgb8_pixel_t(red, green, blue);
		}
	}
	catch (exception& e)
	{
		cout << "Could not parse colour 2. Please provide colour in valid format, for example --colour2 255 255 255";
		return false;
	}
	return true;
}

int main(int argc, char* argv[])
{
	bool shouldContinue = GetParameters(argc, argv);
	if (!shouldContinue)
	{
		return 1;
	}
	shouldContinue = ValidateParameters();
	if (!shouldContinue)
	{
		return 1;
	}
	
	rgb8_image_t img;
	jpeg_read_image(inputImage, img);
	int imageWidth = img.width();
	int imageHeight = img.height();
	printf("Input image width: %i\n", imageWidth);
	printf("Input image height: %i\n", imageHeight);
	vector<vector<quaternion<double>>> quaternionMatrix = ConvertRGBImageToQuaternionMatrix(img);
	vector<vector<quaternion<double>>> quaternionMatrixAfterFilter = ApplyHypercomplexFilter(quaternionMatrix, colour1, colour2);
	if (showVectorPart)
	{
		stringstream ss;
		ss << outputDir << basename(inputImage) << "_VectorPart.jpg";
		string outputImagePath = ss.str();
		rgb8_image_t vectorImage = ConvertQuaternionMatrixToRGBImage(quaternionMatrixAfterFilter);
		try
		{
			jpeg_write_view(outputImagePath, view(vectorImage));
		}
		catch (exception& e)
		{
			cout << "Could not save vector part of quaternion matrix to " << outputImagePath;
			return 1;
		}
	}
	if (showScalarPart)
	{
		stringstream ss;
		ss << outputDir << basename(inputImage) << "_ScalarPart.jpg";
		string outputImagePath = ss.str();
		gray8_image_t scalarImage = ConvertQuaternionMatrixToGrayImage(quaternionMatrixAfterFilter);
		try
		{
			jpeg_write_view(outputImagePath, view(scalarImage));
		}
		catch (exception& e)
		{
			cout << "Could not save scalar part of quaternion matrix to " << outputImagePath;
			return 1;
		}
	}
	stringstream ss;
	ss << outputDir << basename(inputImage) << "_Edges.jpg";
	string outputImagePath = ss.str();
	gray8_image_t edgeImage = ApplyRelaxedThreshold(quaternionMatrixAfterFilter, colour1, colour2);
	try
	{
		jpeg_write_view(outputImagePath, view(edgeImage));
	}
	catch (exception& e)
	{
		cout << "Could not save image containing edges from input image to " << outputImagePath;
		return 1;
	}
	cout << "Done\n";
#ifdef DEBUG
	getchar();
#endif
	return 0;
}


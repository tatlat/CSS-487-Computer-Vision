// Simple program to use the CSS 487 Image class to create a new image that is 
// a scale, shear, rotate, and translate of "test.gif". Uses matrix 
// multiplication and bilinear interpolation over the 4 surrounding pixels to 
// accomplish the transformations. 
// Author: Tanvir Tatla

#define _USE_MATH_DEFINES 
#include <cmath>
#include "Image.h"
#include <iostream>

// coordinate represents the x and y (Cartesian) coordinate in an image/grid
struct coordinate {
	double x;
	double y;
};

// clamp
// Preconditions: N/A
// Postconditions: returns an int that is between 0 and the bound - 1. 
//                 returns n if those conditions are already met, 0 if n is 
//                 negative, or bound - 1 if n is >= bound.
int clamp(const double& n, const int& bound) {
	int m = int(max(0.0, n));
	return min(bound - 1, m);
}

// getNeighborCoordinates
// Preconditions: N/A
// Postconditions: finds the coordinates of the 4 surrounding pixels around 
//                 xPrime and yPrime. assigns those values to tL, tR, bL, and 
//                 bR's x and y fields.
// Parameters: tL (top left coordinate), tR (top right), bL (bottom left), and 
//             bR (bottom right).
void getNeighborCoordinates(coordinate& tL, coordinate& tR, coordinate& bL, 
	coordinate& bR, const double& xPrime, const double& yPrime, int width, 
	int height) {

	tL.x = clamp(xPrime, width), tL.y = clamp(yPrime, height);
	tR.x = ceil(xPrime), tR.y = tL.y;
	tR.x = clamp(tR.x, width);
	bL.x = tL.x, bL.y = ceil(yPrime);
	bL.y = clamp(bL.y, height);
	bR.x = tR.x, bR.y = bL.y;
}

// getNeighborPixels
// Preconditions: input image passed in is correctly allocated/formatted
// Postconditions: assigns the 4 surrounding pixels to tL, tR, bL, and bR 
//                 based on the coordinates passed into the function.
// Parameters: tL (top left pixel), tR (top right), bL (bottom left), 
//             bR (bottom right). tLC (top left coordinate),
void getNeighborPixels(const Image& input, pixel& tL, pixel& tR, pixel& bL, 
	pixel& bR, const coordinate& tLC, const coordinate& tRC, 
	const coordinate& bLC, const coordinate& bRC) {

	tL = input.getPixel(int(tLC.y), int(tLC.x));
	tR = input.getPixel(int(tRC.y), int(tRC.x));
	bL = input.getPixel(int(bLC.y), int(bLC.x));
	bR = input.getPixel(int(bRC.y), int(bRC.x));
}

// calculateWeights
// Preconditions: dX and dY are not zero
// Postconditions: calculates and assigns the weights required for bilinear 
//                 interpolation to w1, w2, w3, and w4.
// Parameters: w1 (weight1), dX (deltaX), dY (deltaY).
void calculateWeights(double& w1, double& w2, double& w3, double& w4, 
	const coordinate& tL, const coordinate& tR, const coordinate& bL, 
	const coordinate& bR, const double& xPrime, const double& yPrime,
	const double& dX, const double& dY) {

	w1 = (tR.x - xPrime) / dX;
	w2 = (xPrime - tL.x) / dX;
	w3 = (bL.y - yPrime) / dY;
	w4 = (yPrime - tL.y) / dY;
}

// getWeightedAverage
// Preconditions: N/A
// Postconditions: returns the weighted average of lhs and rhs.
// Parameters: lhs (value corresponding to weight1), rhs (value corresponding 
//             to w2).
byte getWeightedAverage(const double& w1, const double& w2, 
	const byte& lhs, const byte& rhs) {
	return static_cast<byte>((w1 * lhs) + (w2 * rhs));
}

// interpolate
// Preconditions: input image passed in is correctly allocated/formatted
// Postconditions: assigns the correct red, green, and blue values after doing 
//                 bilinear interpolation
void interpolate(byte& red, byte& green, byte& blue, const double& xPrime, 
	const double& yPrime, const int& width, const int& height, 
	const Image& input) {

	// coordinates for the 4 pixels that surround xPrime, yPrime.
	coordinate topL, topR, bottomL, bottomR; 
	pixel topLeft, topRight, bottomLeft, bottomRight; // neighboring pixels
	byte topRed, bottomRed, topGreen, bottomGreen, topBlue, bottomBlue;
	double deltaX, deltaY, weight1 = 0, weight2 = 0, weight3 = 0, weight4 = 0;

	// find the coordinates of the pixels surrounding (xPrime, yPrime)
	getNeighborCoordinates(topL, topR, bottomL, bottomR, xPrime, yPrime, width,
		height);

	// retrieve the pixels at the coordinates
	getNeighborPixels(input, topLeft, topRight, bottomLeft, bottomRight, topL,
		topR, bottomL, bottomR);

	deltaX = topR.x - topL.x;
	deltaY = bottomL.y - topL.y;

	// avoid divide by zero
	if (deltaX == 0 || deltaY == 0) {
		red = input.getPixel(int(yPrime), int(xPrime)).red;
		green = input.getPixel(int(yPrime), int(xPrime)).green;
		blue = input.getPixel(int(yPrime), int(xPrime)).blue;
		return;
	}

	// calculate the weights neccessary for bilinear interpolation
	calculateWeights(weight1, weight2, weight3, weight4, topL, topR, bottomL, 
		bottomR, xPrime, yPrime, deltaX, deltaY);

	// linear interpolation between the 2 top values
	topRed = getWeightedAverage(weight1, weight2, topLeft.red, topRight.red);
	topGreen = getWeightedAverage(weight1, weight2, topLeft.green, topRight.green);
	topBlue = getWeightedAverage(weight1, weight2, topLeft.blue, topRight.blue);

	// interpolation between 2 bottom values.
	bottomRed = getWeightedAverage(weight1, weight2, bottomLeft.red, 
		bottomRight.red);
	bottomGreen = getWeightedAverage(weight1, weight2, bottomLeft.green, 
		bottomRight.green);
	bottomBlue = getWeightedAverage(weight1, weight2, bottomLeft.blue, 
		bottomRight.blue);

	// interpolation between combined top and bottom values to get final values.
	red = getWeightedAverage(weight3, weight4, topRed, bottomRed);
	green = getWeightedAverage(weight3, weight4, topGreen, bottomGreen);
	blue = getWeightedAverage(weight3, weight4, topBlue, bottomBlue);

	// make sure values are between 0-255 and cast to byte
	red = static_cast<byte>(clamp(red, 255));
	green = static_cast<byte>(clamp(green, 255));
	blue = static_cast<byte>(clamp(blue, 255));
}

// recenterPixel
// Preconditions: N/A
// Postconditions: recenters xPrime and yPrime around the center of the image
void recenterPixel(double& xPrime, double& yPrime, const double& xOrigin, 
	const double& yOrigin) {
	xPrime += xOrigin;
	yPrime += yOrigin;
}

// initializeMatrix
// Preconditions: N/A
// Postconditions: assigns 0 to each row and column in a 2x2 matrix
void initializeMatrix(double *matrix[2]) {
	// loop over each cell in matrix
	for (int row = 0; row < 2; row++) {
		matrix[row] = new double[2]; // add new row

		// initialize each column in row
		for (int col = 0; col < 2; col++)
			matrix[row][col] = 0.0;
	}
}

// deleteMatrix
// Preconditions: matrix is 2x2 array
// Postconditions: frees memory from dynamically allocated 2d arrays
void deleteMatrix(double** matrix) {
	// loop over rows
	for (int i = 0; i < 2; i++)
		delete[] matrix[i]; // delete row

	delete[] matrix;
}

// fillMatrix
// Preconditions: matrix is 2x2
// Postconditions: fills in the matrix using the values of the parameters
void fillMatrix(double** matrix, const double& tLeft, const double& tRight, 
	const double& bLeft, 
	double bRight) {

	matrix[0][0] = tLeft;
	matrix[0][1] = tRight;
	matrix[1][0] = bLeft;
	matrix[1][1] = bRight;
}

// scaleMatrix
// Preconditions: N/A
// Postconditions: returns a 2x2 matrix representing the inverse scale 
//                 transformation. If xScale or yScale are 0, then they
//                 are reassigned a value of 1.
double** scaleMatrix(double& xScale, double& yScale) {
	double** scale = new double* [2];
	initializeMatrix(scale);

	// avoid dividing by zero
	if (xScale == 0) xScale = 1;
	if (yScale == 0) yScale = 1;

	// determinant is 1 / (a*d - b*c)
	double determinant = 1 / (xScale * yScale); // b and c are both zero
	double a = determinant * yScale, d = determinant * xScale;

	fillMatrix(scale, a, 0, 0, d);
	return scale;
}

// shearMatrix
// Preconditions: N/A
// Postconditions: returns a 2x2 matrix representing the inverse shear 
//                 transformation
// Paramters: k is the amount of shearing that should be done to the image.
double** shearMatrix(const double& k) {
	double** shear = new double* [2];
	initializeMatrix(shear);

	// determinant is 1 so no need to muliply
	fillMatrix(shear, 1, -k, 0, 1);
	return shear;
}

// rotationMatrix
// Preconditions: N/A
// Postconditions: returns a 2x2 matrix representing the inverse rotation 
//                 transformation
double** rotationMatrix(const double& angle) {
	double** rotation = new double* [2];
	initializeMatrix(rotation);

	double theta = angle * M_PI / 180; // convert degrees to radians
	double cosine = cos(theta), sine = sin(theta);

	// determinant is 1 so no need to multiply
	fillMatrix(rotation, cosine, sine, -sine, cosine);
	return rotation;
}

// multiplyMatrices
// Preconditions: lhs and rhs are both 2x2 matrices
// Postconditions: returns a 2x2 matrix representing the product of lhs and rhs
double** multiplyMatrices(double** lhs, double** rhs) {
	double** product = new double* [2];
	initializeMatrix(product);
	double tL, tR, bL, bR;

	// calculate the new values that go into our product matrix
	tL = (lhs[0][0] * rhs[0][0]) + (lhs[0][1] * rhs[1][0]);
	tR = (lhs[0][0] * rhs[0][1]) + (lhs[0][1] * rhs[1][1]);
	bL = (lhs[1][0] * rhs[0][0]) + (lhs[1][1] * rhs[1][0]);
	bR = (lhs[1][0] * rhs[0][1]) + (lhs[1][1] * rhs[1][1]);

	fillMatrix(product, tL, tR, bL, bR);
	deleteMatrix(lhs); // destroy since we no longer need lhs and rhs after this
	deleteMatrix(rhs);
	return product;
}

// applyTransformation
// Preconditions: matrix is 2x2 and represents a (combination) transformation
// Postconditions: calculates and assigns new values to xPrime and yPrime 
//                 using matrix multiplication. 
// Parameters: qX and qY are the original x and y (column and row) values to 
//             transform.
void applyTransformation(double** matrix, double& xPrime, double& yPrime, 
	const double& qX, const double& qY, const double& xOrigin, 
	const double& yOrigin) {

	// convert to Cartesian coordinates
	double x = qX - xOrigin, y = qY - yOrigin; 

	xPrime = (x * matrix[0][0]) + (y * matrix[0][1]);
	yPrime = (x * matrix[1][0]) + (y * matrix[1][1]);
}

// transformImage
// Preconditions:  input image passed in is correctly allocated/formatted
// Postconditions: returns a new image that is scaled, sheared, rotated, 
//                 and translated according to the parameters
Image transformImage(const Image& input, double& xScale, double& yScale, 
	const double& angle, const double& shearScale, const double& xTranslate, 
	const double& yTranslate) {

	int height = input.getRows(), width = input.getCols();
	// create black image with same dimensions as input
	Image output(height, width);
	double xOrigin = width / 2.0, yOrigin = height / 2.0; // center coordinates
	byte red, blue, green;

	// get the inverse transformation matrices
	double** scale = scaleMatrix(xScale, yScale);
	double** shear = shearMatrix(shearScale);
	double** rotation = rotationMatrix(angle);

	// combine (multiply) each transformation
	double** partialProduct = multiplyMatrices(scale, shear);
	double** finalMatrix = multiplyMatrices(rotation, partialProduct);

	double xPrime = 0, yPrime = 0, qX = 0, qY = 0;

	// loop over each pixel in the input image
	for (int row = 0; row < height; row++) {
		qY = double(row) - yTranslate; // apply translation first

		for (int col = 0; col < width; col++) {
			qX = double(col) - xTranslate;

			// apply the rest of the transformations
			applyTransformation(finalMatrix, xPrime, yPrime, qX, qY, xOrigin, yOrigin);
			// recenter the pixel so it can fit into the image grid
			recenterPixel(xPrime, yPrime, xOrigin, yOrigin);

			// if xPrime or yPrime go beyond boundaries of the image, skip iteration
			if (xPrime < 0 || xPrime >= width || yPrime < 0 || yPrime >= height) 
				continue;

			// use bilinear interpolation to figure out correct color of pixel
			interpolate(red, green, blue, xPrime, yPrime, width, height, input);
			output.setPixel(row, col, red, green, blue); // set pixel color
		}
	}

	deleteMatrix(finalMatrix);
	return output;
}

// main method
// Preconditions:  test1.gif exists and is a correctly formatted GIF image, 
//                 all arguments (except the name of the program) can be 
//                 converted to valid double values or otherwise the initial 
//                 values will be used. Arguments should in the following 
//                 order: xScale yScale xTransalte yTranslate angle shear
// Postconditions: Creates an image output.gif that is the scale, shear, 
//                 rotation, and translation of test.gif
int main(int argc, char* argv[])
{
	// check if enough arguments were passed
	if (argc < 7) {
		std::cout << "too few arguments" << endl;
		return EXIT_FAILURE;
	}

	// convert and assign arguments to their corresponding double variables
	double xScale = 1.0, yScale = 1.0, xTranslate = 0.0, yTranslate = 0.0, 
		angle = 0.0, shear = 0.0;
	sscanf_s(argv[1], "%lf", &xScale);
	sscanf_s(argv[2], "%lf", &yScale);
	sscanf_s(argv[3], "%lf", &xTranslate);
	sscanf_s(argv[4], "%lf", &yTranslate);
	sscanf_s(argv[5], "%lf", &angle);
	sscanf_s(argv[6], "%lf", &shear);

	// open the input image
	Image input("test1.gif"); 

	// transform the input image according to the arguments
	Image output = transformImage(input, xScale, yScale, angle, shear, 
		xTranslate, yTranslate);

	// write transformed image to output
	output.writeImage("output.gif");
	return 0;
}
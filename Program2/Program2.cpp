// Simple program that uses the 487 Image class to create 2 new images that are
// a smoothing and edge detection of "test2.gif". Uses linear 
// filtering/convolution to accomplish the transformations. Edge detection 
// also uses bilinear interpolation. The smoothing and edge detection kernels 
// are hard-coded but the filter function can work with any kernel.
// Author: Tanvir Tatla

#define _USE_MATH_DEFINES 
#include <cmath>
#include "Image.h"
#include <iostream>

// Sy and Sx are the smoothing kernels
// Ey and Ex are the gradient kernels for edge detection
Image Sy(3, 1), Sx(1, 3), Ey(3, 1), Ex(1, 3);

// Coordinate represents the x and y locations of a pixel in an Image
struct Coordinate {
	int x = 0;
	int y = 0;

	Coordinate() {}

	~Coordinate() {}

	// Constructor
	Coordinate(int nX, int nY) {
		x = nX;
		y = nY;
	}
};

// setKernelX
// Preconditions: k must have 1 row and 3 columns
// Postconditions: sets each cell in the Image kernel
// Parameters: k is the kernel. a, b, c are the values to fill the kernel
//															from left to right.
void setKernelX(Image& k, float a, float b, float c) {
	k.setFloat(0, 0, a);
	k.setFloat(0, 1, b);
	k.setFloat(0, 2, c);
}

// setKernelX
// Preconditions: k must have 3 rows and 1 column
// Postconditions: sets each cell in the Image kernel
// Parameters: k is the kernel. a, b, c are the values to fill the kernel
//															from top to bottom.
void setKernelY(Image& k, float a, float b, float c) {
	k.setFloat(0, 0, a);
	k.setFloat(1, 0, b);
	k.setFloat(2, 0, c);
}

// clamp
// Preconditions: N/A
// Postconditions: returns n if n is between 0 and border. Otherwise returns 0
//									if n < 0, or border - 1 if n >= border
// Parameters: k is the kernel. a, b, c are the values to fill the kernel
//															from left to right.
int clamp(const int& n, int border) {
	if (n < 0) return 0;
	else if (n >= border) return border - 1;
	else return n;
}

// floatingPoint
// Preconditions: input must be a grey scale image
// Postconditions: converts a grey scale image to floating point and returns it
Image floatingPoint(const Image& input) {
	int height = input.getRows(), width = input.getCols();
	Image floating(height, width);

	// for each pixel in the input
	for (int r = 0; r < height; r++) {
		for (int c = 0; c < width; c++) {

			// cast the grey value to float and put the result into floating
			floating.setFloat(r, c, input.getPixel(r, c).grey);
		}
	}
	return floating;
}

// greyScale
// Preconditions: input must be a floating point image
// Postconditions: converts a floating point image to grey scale and returns it
Image greyScale(const Image& input) {
	int height = input.getRows(), width = input.getCols();
	Image grey(height, width);

	// for each pixel in the input
	for (int r = 0; r < height; r++) {
		for (int c = 0; c < width; c++) {

			// cast the float to byte and set the grey value of the output
			grey.setGrey(r, c, static_cast<byte>(input.getFloat(r, c)));
		}
	}
	return grey;
}

// filter
// Preconditions: input must be a floating point image
// Postconditions: convolves input image with the kernel and returns the result
Image filter(const Image& input, const Image& kernel, const int& xCenter, const int& yCenter) {
	int width = input.getCols(), height = input.getRows();
	int kWidth = kernel.getCols(), kHeight = kernel.getRows();
	Image output(height, width);
	int x = 0, y = 0; 
	float sum = 0.0f;

	// for each pixel in the input
	for (int row = 0; row < height; row++) {
		for (int col = 0; col < width; col++) {

			// for each cell in the kernel from last to first
			for (int krow = kHeight - 1; krow >= 0; krow--) {
				for (int kcol = kWidth - 1; kcol >= 0; kcol--) {

					// get x and y of input image
					x = col + (xCenter - kcol);
					y = row + (yCenter - krow);

					// make sure the x and y don't go beyond the borders of the image
					x = clamp(x, width);
					y = clamp(y, height);

					// convolve input pixel with kernel
					sum += input.getFloat(y, x) * kernel.getFloat(krow, kcol);
				}
			}

			// add the sum to the output and reset sum for next pixel
			output.setFloat(row, col, sum);
			sum = 0.0f;
		}
	}

	return output;
}

// smoothImage
// Preconditions: input must be a floating point Image
// Postconditions: smooths the input by convolving with Sx and Sy by the 
//									number of repetitions (parameter)
Image smoothImage(const Image& input, const int& repetitions) {
	Image output = input;

	// Smooth in the x direction repetition times
	for (int i = 0; i < repetitions; i++) {
		output = filter(output, Sx, 1, 0);
	}


	// then again in the y direction
	for (int i = 0; i < repetitions; i++) {
		output = filter(output, Sy, 0, 1);
	}

	return output;
}

// getNeighborCoordinates
// Preconditions: N/A
// Postconditions: finds the coordinates of the 4 surrounding pixels around 
//                 x and y. assigns those values to tL, tR, bL, and bR's 
//                 x and y fields.
// Parameters: tL (top left coordinate), tR (top right), bL (bottom left), and 
//             bR (bottom right).
void getNeighborCoordinates(Coordinate& tL, Coordinate& tR, Coordinate& bL,
	Coordinate& bR, const float& x, const float& y, int width, int height) {

	tL.x = int(clamp(int(x), width)), tL.y = int(clamp(int(y), height));
	tR.x = int(ceil(x)), tR.y = tL.y;
	tR.x = (clamp(tR.x, width));
	bL.x = tL.x, bL.y = int(ceil(y));
	bL.y = int(clamp(bL.y, height));
	bR.x = tR.x, bR.y = bL.y;
}

// getNeighborPixels
// Preconditions: input image passed in is correctly allocated/formatted
// Postconditions: assigns the 4 surrounding pixels to tL, tR, bL, and bR 
//                 based on the coordinates passed into the function.
// Parameters: tL (top left pixel), tR (top right), bL (bottom left), 
//             bR (bottom right). tLC (top left coordinate),
void getNeighborPixels(const Image& input, pixel& tL, pixel& tR, pixel& bL,
	pixel& bR, const Coordinate& tLC, const Coordinate& tRC,
	const Coordinate& bLC, const Coordinate& bRC) {

	tL = input.getPixel(tLC.y, tLC.x);
	tR = input.getPixel(tRC.y, tRC.x);
	bL = input.getPixel(bLC.y, bLC.x);
	bR = input.getPixel(bRC.y, bRC.x);
}

// calculateWeights
// Preconditions: N/A
// Postconditions: calculates and assigns the weights required for bilinear 
//                 interpolation to w1, w2, w3, and w4.
// Parameters: w1 (weight1), dX (deltaX), dY (deltaY).
void calculateWeights(float& w1, float& w2, float& w3, float& w4,
	const Coordinate& tL, const Coordinate& tR, const Coordinate& bL,
	const Coordinate& bR, const float& x, const float& y, const float& dX, 
	const float& dY) {
	
	// avoid divide by zero
	if (dX != 0) {
		w1 = (tR.x - x) / dX;
		w2 = (x - tL.x) / dX;
	}
	else { 
		w1 = 1;
		w2 = 0;
	}

	// avoid divide by zero again
	if (dY != 0) {
		w3 = (bL.y - y) / dY;
		w4 = (y - tL.y) / dY;
	}
	else {
		w3 = 1;
		w4 = 0;
	}
}

// getWeightedAverage
// Preconditions: N/A
// Postconditions: returns the weighted average of lhs and rhs.
// Parameters: lhs (value corresponding to weight1), rhs (value corresponding 
//             to w2).
float getWeightedAverage(const float& w1, const float& w2,
	const float& lhs, const float& rhs) {
	return (w1 * lhs) + (w2 * rhs);
}

// interpolate
// Preconditions: input must be a floating point image
// Postconditions: performs bilinear interpolatation to find the floatVal of 
//                 x and y by using the 4 nearest pixels 
float interpolate(const Image& input, const float& x, const float& y) {
	int height = input.getRows(), width = input.getCols();
	// coordinates for the 4 pixels that surround x, y.
	Coordinate topL, topR, bottomL, bottomR;
	pixel topLeft, topRight, bottomLeft, bottomRight; // 4 neighboring pixels
	float dX, dY, w1, w2, w3, w4; // delta and weights
	float top, bottom, output;

	// find the coordinates of the pixels surrounding (x, y)
	getNeighborCoordinates(topL, topR, bottomL, bottomR, x, y, width, height);

	// retrieve the pixels at the coordinates
	getNeighborPixels(input, topLeft, topRight, bottomLeft, bottomRight, topL,
		topR, bottomL, bottomR);

	// calculate deltas
	dX = float(topR.x - topL.x);
	dY = float(bottomL.y - topL.y);

	// calculate the weights neccessary for bilinear interpolation
	calculateWeights(w1, w2, w3, w4, topL, topR, bottomL,
		bottomR, x, y, dX, dY);

	// get averages for top and bottom half
	top = getWeightedAverage(w1, w2, topLeft.floatVal, topRight.floatVal);
	bottom = getWeightedAverage(w1, w2, bottomLeft.floatVal, bottomRight.floatVal);

	// average top and bottom to get final value
	output = getWeightedAverage(w3, w4, top, bottom);

	// output should be between 0 and 255
	if (output < 0) output = 0.0;
	else if (output > 255.0) output = 255.0;

	return output;
}

// nsMaxima
// Preconditions: all Images must be floating point images and have the exact 
//                same size. r and c must not go beyond the number of rows and 
//                columns in any image. The pixel at r,c in gMag must be an edge 
//                candidate (have floatVal >= 10).
// Postconditions: performs non-maxima suppression to find the edges in gMag
// Parameters: gMag - the gradient magnitude of gradX and gradY. 
//						 gradX - result of convolving an image with Ex (x gradient kernel)
//						 gradY - result of convolving image with Ey (y gradient kernel)
float nsMaxima(const Image& gMag, const Image& gradX, const Image& gradY, int r, int c) {
	float gX = gradX.getFloat(r, c) / gMag.getFloat(r, c);
	float gY = gradY.getFloat(r, c) / gMag.getFloat(r, c);

	// calculate the coordinates of points R and P
	float rX = c + gX, rY = r + gY, pX = c - gX, pY = r - gY; 

	// interpolate to get the floatVal of points R and P
	float gR = interpolate(gMag, rX, rY);
	float gP = interpolate(gMag, pX, pY);

	float gQ = gMag.getFloat(r, c);

	// check if r,c in gMag is an edge
	if (gQ > gR && gQ > gP) return 255.0;
	return 0.0;
}

// getMagnitude
// Preconditions: x must be an x-direction gradient, y must be y-direction gradient
//								of the same image as x
// Postconditions: returns the magnitude gradient image of x and y 
Image getMagnitude(const Image& x, const Image& y) {
	int height = x.getRows(), width = x.getCols();
	Image mag(height, width);
	float xSquare, ySquare, root;

	// for each pixel in the gradients
	for (int r = 0; r < height; r++) {
		for (int c = 0; c < width; c++) {

			// calculate the magnitude
			xSquare = pow(x.getFloat(r, c), 2);
			ySquare = pow(y.getFloat(r, c), 2);
			root = sqrt(xSquare + ySquare);
			mag.setFloat(r, c, root); // place into magnitude gradient
		}
	}

	return mag;
}

// edgeDetection
// Preconditions: input must be a floating point image
// Postconditions: returns a floating point image with the edges highlighted
//									in white and everything else is black.
Image edgeDetection(const Image& input) {
	int height = input.getRows(), width = input.getCols();
	Image output(height, width);
	Image gradientX, gradientY, gMag;

	gradientX = filter(input, Ex, 1, 0); // get the x-direction gradient
	gradientY = filter(input, Ey, 0, 1); // same for y-direction
	gMag = getMagnitude(gradientX, gradientY); // get magnitude gradient

	// for each pixel in the magnitude gradient
	for (int row = 0; row < height; row++) {
		for (int col = 0; col < width; col++) {

			// if the pixel is an edge candidate
			if (gMag.getFloat(row, col) >= 10.0) {

				// perform non-maxima suppression
				float q = nsMaxima(gMag, gradientX, gradientY, row, col);
				output.setFloat(row, col, q);
			}
			else {
				output.setFloat(row, col, 0.0); // otherwise make the pixel black
			}
		}
	}

	return output;
}

// main method
// Preconditions:  test2.gif exists and is a correctly formatted grey scale 
//                 GIF image, the argument can be read as an int value. 
// Postconditions: Creates 2 images: smooth.gif and edges.gif. The former 
//									is the result of smoothing test2.gif with linear filtering.
//									The former is the result of edge detection also using convolution.
int main(int argc, char* argv[])
{
	// check if enough arguments were passed
	if (argc < 2) {
		std::cout << "too few arguments" << endl;
		return EXIT_FAILURE;
	}

	// read number of repetitions from argument list
	int repetitions = 0;
	sscanf_s(argv[1], "%d", &repetitions);

	// create the smoothing and edge-detection gradient kernels
	setKernelX(Sx, 0.25, 0.5, 0.25);
	setKernelY(Sy, 0.25, 0.5, 0.25);
	setKernelX(Ex, -1.0, 0.0, 1.0);
	setKernelY(Ey, -1.0, 0.0, 1.0);

	// open the input image and convert to floating point
	Image input("test2.gif");
	Image fInput = floatingPoint(input);

	// perform smoothing and edge detection
	Image smooth = smoothImage(fInput, repetitions);
	Image edges = edgeDetection(smooth);

	// convert back to greyScale
	Image gSmooth = greyScale(smooth);
	Image gEdges = greyScale(edges);

	// write output to disk
	gSmooth.writeGreyImage("smooth.gif");
	gEdges.writeGreyImage("edges.gif");
	return 0;
}
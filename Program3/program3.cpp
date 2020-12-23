// program3.cpp
// This program modifies an image named "foreground.jpg" by replacing the most 
// common color with pixels from "background.jpg". The effect is similar to 
// green screen or chroma key. This code also modifies background.jpg by 
// flipping it horizontally, converting it to greyscale, blurring it, and 
// detecting edges using OpenCV methods. Furthermore, it modifies test3.jpg by 
// warping the perspective, applying a color map, and adding text over the image.
// Author: Tanvir Tatla

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

// color holds the rgb values for a pixel
struct color {
	int red = 0;
	int green = 0;
	int blue = 0;
};

// createHistogram - traverses the pixels in foreground and tallies the votes
//		for colors in their respective buckets.
// preconditions: foreground is a color byte image, size is greater than zero
// postconditions: a 3d integer matrix that represents a color histogram
Mat createHistogram(const Mat& foreground, const int& size, const int& bucketSize) {

	// create an array of the histogram dimensions
	// size is a constant - the # of buckets in each dimension
	int dims[] = { size, size, size };
	// create 3D histogram of integers initialized to zero	
	Mat hist(3, dims, CV_32S, Scalar::all(0));

	// loop over pixels in foreground
	for (int row = 0; row < foreground.rows - 1; row++) {
		for (int col = 0; col < foreground.cols - 1; col++) {

			// get the pixel at current index
			Vec3b pixel = foreground.at<Vec3b>(row, col);
			uchar blue = pixel.val[0]; // get rgb values of pixel
			uchar green = pixel.val[1];
			uchar red = pixel.val[2];

			// determine which bucket to increment
			int r = red / bucketSize; 
			int g = green / bucketSize;
			int b = blue / bucketSize;

			// increment bucket
			hist.at<int>(r, g, b)++;
		}
	}

	return hist;
}

// getCommonColor - traverses a color histogram
// preconditions: color histogram is 3d integer matrix, size is greater than 0
// postconditions: returns color containing the rgb combination with the most 
//		votes
color getCommonColor(const Mat& hist, const int& size, const int& bucketSize) {
	color output;
	int max = 0; // max number of votes

	// loop over color histogram
	for (int r = 0; r < size; r++) {
		for (int g = 0; g < size; g++) {
			for (int b = 0; b < size; b++) {

				// get number of votes at current index
				int temp = hist.at<int>(r, g, b);

				// if votes greater than max
				if (temp > max) {
					max = temp; // set max to new value
					output.red = r; // record bin index
					output.green = g;
					output.blue = b;
				}
			}
		}
	}

	// convert bin index to most common color
	output.red = output.red * bucketSize + bucketSize / 2;
	output.green = output.green * bucketSize + bucketSize / 2;
	output.blue = output.blue * bucketSize + bucketSize / 2;
	
	return output;
}

// closeToColor - Used to check if the color of the current pixel in foreground
//		is no more than bucketSize away from the most common color
// preconditions: N/A
// postconditions: returns true if color curr is within bucketSize away from 
//		color c, returns false if at least one color band is more than
//		bucketSize away
bool closeToColor(const color& c, const color& curr, const int& bucketSize) {
	// check each color band individually to see if it is more than bucketSize 
	// away
	if (curr.red > c.red + bucketSize || curr.red < c.red - bucketSize) 
		return false;
	if (curr.green > c.green + bucketSize || curr.green < c.green - bucketSize) 
		return false;
	if (curr.blue > c.blue + bucketSize || curr.blue < c.blue - bucketSize) 
		return false;
	return true;
}

// replaceColor - replaces the pixels in foreground that match with the most 
//		common color of color. The pixels from background are used 
// preconditions: input is a color byte image, size is greater than 0
// postconditions: a new image with the same size of input but some pixels are 
//		replaced with pixels of background
Mat replaceColor(const Mat& input, const Mat& background, const color& common, 
	const int& bucketSize) {

	Mat output = input.clone();
	int x = 0, y = 0;

	// loop over pixels of foreground
	for (int r = 0; r < output.rows - 1; r++) {
		for (int c = 0; c < output.cols - 1; c++) {

			Vec3b pixel = output.at<Vec3b>(r, c); // get current pixel
			color curr;
			curr.blue = pixel.val[0]; // get color of current pixel
			curr.green = pixel.val[1];
			curr.red = pixel.val[2];

			// check if current pixel's color is close to common
			if (closeToColor(common, curr, bucketSize)) {
				// replace foreground pixel with corresponding background pixel
				y = r % background.rows;
				x = c % background.cols;
				output.at<Vec3b>(r, c) = background.at<Vec3b>(y, x);
			}
		}
	}

	return output;
}

// transformImage - flips, converts to greyscale, blurs, and detects edges of 
//		input
// preconditions: input is a color byte image
// postconditions: a new image with the same size is returned after doing all
//		4 operations
Mat transformImage(const Mat& input) {
	Size s = Size(7, 7); // size for Gaussian Blue
	double sX = 2.0, sY = 2.0; // sigma X and Y used for Gaussian Blur
	double thresh1 = 20, thresh2 = 60; // thresholds for canny edge detection
	Mat output = input.clone();

	// Performs operations using OpenCV methods
	flip(output.clone(), output, 1);
	cvtColor(output.clone(), output, COLOR_BGR2GRAY);
	GaussianBlur(output.clone(), output, s, sX, sY);
	Canny(output.clone(), output, thresh1, thresh2);

	return output;
}

// transformImage - warps the perspective, applys a deep green color map, 
//		and places text over the input image
// preconditions: input is grey scale or color byte image
// postconditions: a new image with the same size is returned after doing all
//		3 operations
Mat myTransformImage(const Mat& input) {
	Mat output = input.clone();
	Point2f src[4], dst[4]; // points of quadrilateral

	// The 4 points of quadilateral selected in the src image
	src[0] = Point2f(-45, -45);
	src[1] = Point2f(output.cols + 25, -25);
	src[2] = Point2f(output.cols + 50, output.rows + 25);
	src[3] = Point2f(-25, output.rows + 25);

	// The 4 points where the mapping should go
	dst[0] = Point2f(0, 0); // top left corner
	dst[1] = Point2f(output.cols - 1, 0); // top right
	dst[2] = Point2f(output.cols - 1, output.rows - 1); // bottom right
	dst[3] = Point2f(0, input.rows - 1); // bottom left

	// Perspective Transform Matrix
	Mat pTM = getPerspectiveTransform(src, dst);
	// apply perspective transformation to image
	warpPerspective(output.clone(), output, pTM, output.size());

	// remap the color of the image to green
	applyColorMap(output.clone(), output, COLORMAP_DEEPGREEN);

	Scalar gold = CV_RGB(255, 215, 0); // color of text
	double scale = 4.0; // font size
	int thickness = 7; // font thickness
	Point org = Point(output.cols / 3, output.rows / 5); // location of font

	// lay string as text over image
	putText(output, "Epic", org, FONT_HERSHEY_DUPLEX, scale, gold, thickness);
	return output;
}

// main - reads foreground and background and writes a new image with a chroma
//		key effect using pixels from background to replace the most common color 
//		foreground. Also flips, converts to greyscale, blurs, and detects edges
//		of background. Warps the perspective, applies a color map, and places 
//		text over test3.jpg.
// precondition: foreground.jpg, background.jpg, and test3.jpg exist in the 
//		code directory and are valid JPGs
// postconditions: The image with the chroma key effect is displayed and written
//		to the disk as overlay.jpg. The flipped, greyscale, blurred, and edge image
//		is displayed on the screen and written to disk as output.jpg. The warped
//		perspective image is displayed and written to disk as myoutput.jpg.
int main(int argc, char* argv[])
{	
	const int size = 4;
	const int bucketSize = 256 / size;

	// open foreground and background images
	Mat foreground = imread("foreground.jpg");	
	Mat background = imread("background.jpg");

	// create color histogram of foreground
	Mat hist = createHistogram(foreground, size, bucketSize);
	// get the most common color from the color histogram
	color common = getCommonColor(hist, size, bucketSize);
	// replace the most common color in foreground using pixels from background
	Mat overlay = replaceColor(foreground, background, common, bucketSize);

	namedWindow("overlay", WINDOW_NORMAL);
	resizeWindow("overlay", overlay.cols / 8, overlay.rows / 8);
	imshow("overlay", overlay); // display the overlay image
	waitKey(0);
	imwrite("overlay.jpg", overlay); // write to disk

	// apply the flip, greyscale, blur, and edge detection to background
	Mat output = transformImage(background);
	imshow("output", output); // display the image
	waitKey(0);
	imwrite("output.jpg", output); // write to disk

	// read the original image
	Mat original = imread("test3.jpg");
	// warp the perspective, change color, and place text
	Mat myoutput = myTransformImage(original);
	imshow("myoutput", myoutput); // display the result
	waitKey(0);
	imwrite("myoutput.jpg", myoutput); // write to disk
	return 0;
}
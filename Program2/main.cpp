// Program 2






#include "Image.h"
#include "math.h"
#include <iostream>

void toFloat(Image& input);
void toGrey(Image& input);
Image convolution(const Image& input, const Image& kernel, int rowOrigin, int colOrigin);
Image create_Gmag(const Image& Gy, const Image& Gx);
Image nonMaxSup(const Image& Gy, const Image& Gx, const Image& Gmag);
float bilinearInterpolation(float alpha, float beta, float topLeft, float topRight, float botLeft, float botRight);


int main()
{
	Image output("test2.gif");
	toFloat(output);




	Image kernel1(1, 3);
	Image kernel2(3, 1);
	Image edgeKernel1(1, 3);
	Image edgeKernel2(3, 1);

	kernel1.setFloat(0, 0, 0.25);
	kernel1.setFloat(0, 1, 0.5);
	kernel1.setFloat(0, 2, 0.25);

	kernel2.setFloat(0, 0, 0.25);
	kernel2.setFloat(1, 0, 0.5);
	kernel2.setFloat(2, 0, 0.25);

	edgeKernel1.setFloat(0, 0, -1.0f);
	edgeKernel1.setFloat(0, 1, 0.0f);
	edgeKernel1.setFloat(0, 2, 1.0f);

	edgeKernel2.setFloat(0, 0, -1.0f);
	edgeKernel2.setFloat(1, 0, 0.0f);
	edgeKernel2.setFloat(2, 0, 1.0f);

	int smoothingAmount = 2;
	for (int i = 0; i < smoothingAmount; i++)
	{
		output = convolution(output, kernel1, 0, 1);
		output = convolution(output, kernel2, 1, 0);
	}

	Image Gx = convolution(output, edgeKernel1, 0, 1);
	Image Gy = convolution(output, edgeKernel2, 1, 0);

	Gy.writeFloatImage("myGY.gif");
	Gx.writeFloatImage("myGX.gif");

	toGrey(output);
	output.writeGreyImage("smoothed.gif");


	Image Gmag = create_Gmag(Gy, Gx);

	Gmag.writeFloatImage("Gmag.gif");
	Image c = nonMaxSup(Gy, Gx, Gmag);
	toGrey(c);
	c.writeGreyImage("d.gif");
	c.writeFloatImage("aasdf.gif");
	exit(0);
	/*
	

	
	

	
	
//	output = convolution(output, kernel1, 0, 1);
//	output = convolution(output, kernel1, 0, 1);
//
//	output = convolution(output, kernel2, 1, 0);
//	output = convolution(output, kernel2, 1, 0);

	Image Gx = convolution(output, edgeKernel1, 0, 1);
	Image Gy = convolution(output, edgeKernel2, 1, 0);
	//output = toGrey(input);
	toFloat(output);
	output.writeFloatImage("smoothed.gif");
	Gy.writeFloatImage("myGY.gif");
	Gx.writeFloatImage("myGX.gif");

	Image Gmag = create_Gmag(Gy, Gx);
	
	Image c = nonMaxSup(Gy, Gx, Gmag);
	
	c.writeFloatImage("c.gif");
	
	
	
	
	// take image
	// turn from grey to float
	// smooth x amount of times
	// edge gradient for x and y
	// compute Gmag with distance formaula check disc if u don't remember the convo
	// compute if edge using the non max suppression ish and the interpolation
	// turn from float to grey

	*/
}


void toFloat(Image& input)
{
	Image output(input.getRows(), input.getCols());

	for (int y = 0; y < input.getRows(); y++)
	{
		for (int x = 0; x < input.getCols(); x++)
		{
			input.setFloat(y, x, input.getPixel(y, x).grey);
		}
	}

	//return output;
}

void toGrey(Image& input)
{
	//Image output(input.getRows(), input.getCols());

	for (int y = 0; y < input.getRows(); y++)
	{
		for (int x = 0; x < input.getCols(); x++)
		{
			input.setGrey(y, x, input.getFloat(y, x));
		}
	}

	//return output;
}


Image convolution(const Image& input, const Image& kernel, int rowOrigin, int colOrigin)
{
	Image output(input.getRows(), input.getCols());

	int kCol = kernel.getCols();
	int kRow = kernel.getRows();

	for (float rowIndex = 0.0; rowIndex <= input.getRows() - 1; rowIndex++)
	{
		for (float colIndex = 0.0; colIndex <= input.getCols() - 1; colIndex++)
		{
			float sum = 0.0f;
			
			for (float kRowIndex = kRow - 1; kRowIndex >= 0; kRowIndex--)
			{
				for (float kColIndex = kCol - 1; kColIndex >= 0; kColIndex--)
				{
					int x = colIndex + colOrigin - kColIndex;
					int y = rowIndex + rowOrigin - kRowIndex;

					if (x < 0)
					{
						x = 0;
					}
					else if (x >= input.getCols())
					{
						x = input.getCols() - 1;
					}
					
					if (y < 0)
					{
						y = 0;
					}
					else if (y >= input.getRows())
					{
						y = input.getRows() - 1;
					}

					sum += kernel.getFloat(kRowIndex, kColIndex) * input.getFloat(y, x);
				}

			}
			output.setFloat(rowIndex, colIndex, sum);
			sum = 0.0;
		}
	}

	return output;
}

Image create_Gmag(const Image& Gy, const Image& Gx)
{
	Image output(Gy.getRows(), Gy.getCols());
	for (int rowIndex = 0; rowIndex < Gy.getRows(); rowIndex++)
	{
		for (int colIndex = 0; colIndex < Gy.getCols(); colIndex++)
		{
			output.setFloat(rowIndex, colIndex, sqrtf((Gy.getFloat(rowIndex, colIndex) * Gy.getFloat(rowIndex, colIndex)) +
				(Gx.getFloat(rowIndex, colIndex) * Gx.getFloat(rowIndex, colIndex))));
		}
	}

	return output;
}

int clamp(int n, int border) {
	if (n < 0) return 0;
	if (n > border) return border - 1;
	return n;
}


Image nonMaxSup(const Image& Gy, const Image& Gx, const Image& Gmag)
{
	Image output(Gy.getRows(), Gy.getCols());

	for (float rowIndex = 0.0f; rowIndex < Gy.getRows(); rowIndex++)
	{
		for (int colIndex = 0.0f; colIndex < Gy.getCols(); colIndex++)
		{
			if (Gmag.getFloat(rowIndex, colIndex) >= 10.0f)
			{
			//	Gmag(r + Gy.floatVal((r + Gy(r, c) / Gmag(r, c)), c + Gx.floatVal(r, (c + Gx(r, c)) / Gmag(r, c)))


			// +1
				float alpha = Gx.getFloat(rowIndex, colIndex) / Gmag.getFloat(rowIndex, colIndex);
			//	float plusOneY = Gy.getFloat(rowIndex + Gy.getFloat(rowIndex, colIndex) / Gmag.getFloat(rowIndex, colIndex), colIndex);

				float beta = Gy.getFloat(rowIndex, colIndex) / Gmag.getFloat(rowIndex, colIndex);

				if (alpha == 0.0f && beta == 0.0f)
				{
					output.setFloat(rowIndex, colIndex, 255.0f);
				}
				else
				{

					//float Gr = 500.0f;
					//float Gp = 500.0f;
					//std::cout << rowIndex << " | " << colIndex << " | " << plusOneX << " | " << plusOneY << std::endl;
					float Gr = 500.0f;
					float Gp = 500.0f;




					int xLeft = clamp(floor(colIndex + alpha), Gmag.getCols());
					int xRight = clamp(ceil(colIndex + alpha), Gmag.getCols());
					int yTop = clamp(floor(rowIndex + beta), Gmag.getRows());
					int yBot = clamp(ceil(rowIndex + beta), Gmag.getRows());

					/*
					float Gr = ((1.0f - alpha) * (1.0f - beta) * Gmag.getFloat(floor(rowIndex + beta), floor(colIndex + alpha)) +
						((alpha) * (1.0f - beta) * Gmag.getFloat(floor(rowIndex + beta), ceil(colIndex + alpha))) +
						((1.0f - alpha) * (beta)*Gmag.getFloat(ceil(rowIndex + beta), floor(colIndex + alpha))) +
						((alpha) * (beta)* Gmag.getFloat(ceil(rowIndex + beta), ceil(colIndex + alpha))));

					float Gp = ((1.0f - alpha) * (1.0f - beta) * Gmag.getFloat(floor(rowIndex - beta), floor(colIndex - alpha)) +
						((alpha) * (1.0f - beta) * Gmag.getFloat(floor(rowIndex - beta), ceil(colIndex - alpha))) +
						((1.0f - alpha) * (beta)*Gmag.getFloat(ceil(rowIndex - beta), floor(colIndex - alpha))) +
						((alpha) * (beta)*Gmag.getFloat(ceil(rowIndex - beta), ceil(colIndex - alpha))));



						*/
						//	int 
					float alphaPos = 0.0f;
					float betaPos = 0.0f;

					if (alpha < 0.0)
					{
						alphaPos = -1.0 * alpha;
					}
					if (beta < 0.0f)
					{
						betaPos = -1.0 * beta;
					}

					Gr = bilinearInterpolation(fabs(alpha), fabs(beta),
						Gmag.getFloat(yTop, xLeft),
						Gmag.getFloat(yTop, xRight),
						Gmag.getFloat(yBot, xLeft),
						Gmag.getFloat(yBot, xRight));


					xLeft = clamp(floor(colIndex - alpha), Gmag.getCols());
					xRight = clamp(ceil(colIndex - alpha), Gmag.getCols());
					yTop = clamp(floor(rowIndex - beta), Gmag.getRows());
					yBot = clamp(ceil(rowIndex - beta), Gmag.getRows());

					Gp = bilinearInterpolation(fabs(beta), fabs(beta),
						Gmag.getFloat(yTop, xLeft),
						Gmag.getFloat(yTop, xRight),
						Gmag.getFloat(yBot, xLeft),
						Gmag.getFloat(yBot, xRight));




					/*

					if (beta == 0.0f && alpha != 0.0f)
					{
						Gr = Gmag.getFloat(rowIndex, colIndex + 1);
						Gp = Gmag.getFloat(rowIndex, colIndex - 1);
					}
					else if (beta != 0.0f && alpha == 0.0f)
					{
						Gr = Gmag.getFloat(rowIndex + 1, colIndex);
						Gp = Gmag.getFloat(rowIndex - 1, colIndex);
					}
					*/

					if (Gmag.getFloat(rowIndex, colIndex) > Gr && Gmag.getFloat(rowIndex, colIndex) > Gp)
					{
						output.setFloat(rowIndex, colIndex, 255.0f);
					}
					else
					{
						output.setFloat(rowIndex, colIndex, 0.0f);
					}

				}
				
			}
			else
			{
				output.setFloat(rowIndex, colIndex, 0.0f);
			}
		}
	}

	return output;
}

float bilinearInterpolation(float alpha, float beta, float topLeft, float topRight, float botLeft, float botRight)
{
	float newColor = ((1.0f - alpha) * (1.0f - beta) * topLeft) +
		((alpha) * (1.0f - beta) * botLeft) +
		((1.0f - alpha) * (beta)*topRight) +
		((alpha) * (beta)*botRight);

	return newColor;
}


/*

float newColor = ((1.0f - alpha) * (1.0f - beta) * topLeft) +
		((alpha) * (1.0f - beta) * botLeft) +
		((1.0f - alpha) * (beta) * topRight) +
		((alpha) * (beta) * botRight);
*/
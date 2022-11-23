#include<stdio.h>
#include<iostream>
#include<opencv2/opencv.hpp>
using namespace std;
using namespace cv;

void resize_bilinear(Mat &original_img, Mat &resized_img, double rate) {
	int x, y;
	int original_x, original_y;

	for (y = 0; y < resized_img.rows; y++) {
		for (x = 0; x < resized_img.cols; x++) {
			original_x = (int)(x / rate);
			original_y = (int)(y / rate);

			double fx1 = (double)x / (double)rate - (double)original_x;
			double fx2 = 1 - fx1;
			double fy1 = (double)y / (double)rate - (double)original_y;
			double fy2 = 1 - fy1;

			Vec3b p1 = original_img.at<Vec3b>(original_y, original_x);
			Vec3b p2 = original_img.at<Vec3b>(original_y, original_x + 1);
			Vec3b p3 = original_img.at<Vec3b>(original_y + 1, original_x);
			Vec3b p4 = original_img.at<Vec3b>(original_y + 1, original_x + 1);

			double w1 = fx2 * fy2;
			double w2 = fx1 * fy2;
			double w3 = fx2 * fy1;
			double w4 = fx1 * fy1;

			resized_img.at<Vec3b>(y, x) = w1 * p1 + w2 * p2 + w3 * p3 + w4 * p4;
		}
	}
}

int main()
{
	Mat img;
	img = imread("dog.jpg"); 
	int w = img.cols;
	int h = img.rows;
	if (img.empty()) {
		cout << "No image\n";
	}

	imshow("Original image", img);

	Mat resized_img(h, w, CV_8UC3, Scalar(0));
	resize_bilinear(img, resized_img, 1.2);
	imshow("Resized image", resized_img);
	waitKey(0);
}
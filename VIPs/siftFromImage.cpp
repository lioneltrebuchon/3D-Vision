#include <opencv2/opencv.hpp>
#include <opencv2/xfeatures2d.hpp>
#include <stdio.h>
#include <stdarg.h>
#include <iostream>
#include <fstream> 
#include <vector>
#include <math.h>  
// #include <opencv2/nonfree/features2d.hpp>

using namespace cv;
using namespace std;

// Need to compile with contrib package, sorry: https://github.com/itseez/opencv_contrib
// Compile with g++ siftFromImage.cpp -g -I/usr/local/include/opencv -I/usr/local/include -L/usr/local/lib -lopencv_stitching -lopencv_objdetect -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_ml -lopencv_imgproc -lopencv_flann -lopencv_core -lopencv_xfeatures2d

void calcSIFTDescriptor(string filename);

int main(int argc, char **argv) {
    if (argc != 2) {
        cout << "Usage: ./a.out imagefilename.jpg \n";
        return 0;
    }
    calcSIFTDescriptor(argv[1]);
    return 1;
}

void calcSIFTDescriptor(string filename)
{
    Mat img = imread(filename, CV_LOAD_IMAGE_COLOR);
    float h = img.rows;
    float w = img.cols;
    float scl = h/10;
    Mat graySquare;
    cvtColor(img, graySquare, CV_RGB2GRAY);

    Point2f pt(cvRound(h/2), cvRound(w/2));
    KeyPoint kp(pt, scl);
    std::vector<KeyPoint> kps;
    kps.push_back(kp);

    Mat descriptors;
    cv::Ptr<Feature2D> extractor = xfeatures2d::SIFT::create();
    extractor->compute(graySquare, kps, descriptors);

    Mat img_kps;

    for (int i = 0; i < 128; i++){
        cout << descriptors.at<float>(i);
        if (i != 127) {
            cout <<",";
        } else {
            cout <<"\n";
        }
    }
    // cout << descriptors << "\n";
    
}
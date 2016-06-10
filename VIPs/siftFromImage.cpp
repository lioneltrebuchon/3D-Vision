#include <opencv2/opencv.hpp>
#include <iostream>
#include <stdarg.h>

using namespace cv;
using namespace std;

void calcSIFTDescriptor( Mat img, Mat corners, float scl)
{

    double minX, minY, maxX, maxY;
    minMaxLoc(corners.row(0),&minX, &maxX);
    minMaxLoc(corners.row(1),&minY, &maxY);
    double subX = 0;
    double subY = 0;
    if (minX < 0) {
        subX = minX;
    }
    if (minY < 0) {
        subY = minY;
    }
    for (int i = 0; i < corners.cols; i++){
        corners.at<double>(0,i) = corners.at<double>(0,i) - subX;
        corners.at<double>(1,i) = corners.at<double>(1,i) - subY;
    }
    Mat sortedCorners = corners.clone();
    cv::sort(corners, sortedCorners, CV_SORT_EVERY_ROW);
    float size = min(sortedCorners.at<double>(0,2) - sortedCorners.at<double>(0,1), 
        sortedCorners.at<double>(1,2) - sortedCorners.at<double>(1,1));
    if (size == 0){
        return;
    }
    Mat square = img(Rect(sortedCorners.at<double>(0,2), sortedCorners.at<double>(1,2), size, size));
    Mat graySquare;
    cvtColor(square, graySquare, CV_RGB2GRAY);

    Point2f pt(cvRound(size/2), cvRound(size/2));
    KeyPoint kp(pt, scl);
    std::vector<KeyPoint> kps;
    kps.push_back(kp);

    Mat descriptors;
    DescriptorExtractor* extractor;
    extractor = new DescriptorExtractor();
    extractor->compute(graySquare, kps, descriptors);

    Mat img_kps;
    namedWindow( "Display window", WINDOW_AUTOSIZE );// Create a window for display.
    imshow( "Display window", graySquare );                   // Show our image inside it.

    waitKey(0);
    drawKeypoints(graySquare,kps,img_kps,Scalar::all(-1), DrawMatchesFlags::DEFAULT );

    namedWindow( "Display window", WINDOW_AUTOSIZE );// Create a window for display.
    imshow( "Display window", img_kps );                   // Show our image inside it.

    waitKey(0);

    cout << descriptors << "\n";
    
}
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fcntl.h>
#include <sys/stat.h> 
#include <fstream> 
#include <cstdlib>                      // C standard library
#include <cstdio>                       // C I/O (for sscanf)
#include <cstring>                      // string manipulation
#include <string>                       // for stoi
#include <ANN/ANN.h>                    // ANN declarations
#include <sstream> 
#include <vector>
#include <math.h>  
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

class SiftFeature
{
public:
    float point[3];
    float size;
    float orientation;
    char description[128];
    int cameraIndex;
};

class Camera
{
   public:
      string name;   
      double focalLength; 
      double quaternion[4]; 
      double center[3];
      double radialDistortion;
      vector<SiftFeature> SiftVector; 
};

class sparseSiftFeature
{
public:
    SiftFeature Sift;  
    double modelXY[2];
    double translation[3];
    double rotation[3][3];
    double H[3][3];
    double R1[3][3];
};

class sparseModelPoint{
public:
    double point[3];
    double normal[3];
    vector<sparseSiftFeature> features;
};

enum
{
    SIFT_NAME= ('S'+ ('I'<<8)+('F'<<16)+('T'<<24)),
    MSER_NAME= ('M'+ ('S'<<8)+('E'<<16)+('R'<<24)),
    RECT_NAME= ('R'+ ('E'<<8)+('C'<<16)+('T'<<24)), 
    SIFT_VERSION_4=('V'+('4'<<8)+('.'<<16)+('0'<<24)),
    SIFT_VERSION_5=('V'+('5'<<8)+('.'<<16)+('0'<<24)),
    SIFT_EOF = (0xff+('E'<<8)+('O'<<16)+('F'<<24)),
};

static inline int IsValidVersionName(int value)
{
    return value == SIFT_VERSION_4 || value == SIFT_VERSION_5;
}

static inline int IsValidFeatureName(int value)
{
    return value == SIFT_NAME || value == MSER_NAME;
}

vector<SiftFeature> readSIFT(string filename, int index);
void computeTranslation(Camera c, sparseSiftFeature *s, sparseModelPoint smp);
void computeRotation(Camera c, sparseSiftFeature *s, sparseModelPoint smp);
void computeHomography(Camera c, sparseSiftFeature *s, double normal[3]);
void createVIP(Camera c, string imageName, sparseSiftFeature *s, string patchName);
Mat MakeRotationMatrix(sparseModelPoint smp);  

/******************************
* I am a file that sets up models/features for creating Viewpoint Invariant Patches
* If you don't compile me with the ann library I will be sad :( (See ann user guide)
* Also needs to be compiled with OpenCV now

Compile command: (Probably don't need all the OpenCV libraries listed here)
g++ warpSIFT.cpp -I/path/to/ann_1.1.2/include -L/path/to/ann_1.1.2/lib -lANN -g -I/usr/local/include/opencv -I/usr/local/include -L/usr/local/lib -lopencv_shape -lopencv_stitching -lopencv_objdetect -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_ml -lopencv_imgproc -lopencv_flann -lopencv_core
Lio: g++ -std=c++11  warpSIFT.cpp -I/home/lionelt/ann_1.1.2/include -L/home/lionelt/ann_1.1.2/lib -lANN -g -I/usr/local/include/opencv -I/usr/local/include -L/usr/local/lib -lopencv_shape -lopencv_stitching -lopencv_objdetect -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_ml -lopencv_imgproc -lopencv_flann -lopencv_core


To run, something like that:
./a.out -df ../ETH_example_3D_model/dense3.nvm.cmvs/00/models/option-0000.ply -sf ../ETH_example_3D_model/dense3.nvm -sift ../ETH_example_3D_model/eth_photos_night

* Must be given 3 flags:
*       -df: path to dense file
*       -sf: path to sparse file
*       -sift: path to folder containing images and SIFT files
*              that generated the model. (Plz no '/' at end of path)
*******************************/
int main(int argc, char **argv)
{
    ifstream denseStream;
    ifstream sparseStream;
    string siftFolder;
    int maxPoints = 1024;
    char line[1024];

    // Read in command line arguments and open the file streams
    int i = 1;
    while (i < argc) {                       
        if (!strcmp(argv[i], "-df")) {    
            denseStream.open(argv[++i], ios::in);
            if (!denseStream) {
                cerr << "Cannot open dense file\n";
                exit(1);
            }             
        }
        else if (!strcmp(argv[i], "-sf")) {   
            sparseStream.open(argv[++i], ios::in);
            if (!sparseStream) {
                cerr << "Cannot open sparse file\n";
                exit(1);
            }
        }
        else if (!strcmp(argv[i],"-sift")){
            siftFolder = argv[++i];
        } 
        else if (!strcmp(argv[i],"-n")){
            maxPoints = stoi(argv[++i]);
        }
        else {                                  
            cerr << "Unrecognized option.\n";
            exit(1);
        }
        i++;
    }

    // Store cameras and sift features
    // Creates a vector of Camera objects. Each Camera object
    // has image name, focal length, quaternion, center, distortion
    // and a vector of all SIFT features of that camera/image
    sparseStream.getline(line, 1028);
    sparseStream.getline(line, 1028);
    sparseStream.getline(line, 1028);
    int numCameras = stoi(line);
    vector<Camera> cameras(numCameras);
    // Read a camera line and store all camera parameters
    for (int i = 0; i < numCameras; i++){
        sparseStream.getline(line, 1028);
        istringstream ss(line);
        Camera c;
        ss >> c.name;
        ss >> c.focalLength;
        for (int j = 0; j < 4; j++){
            ss >> c.quaternion[j];
        }
        for (int j = 0; j < 3; j++){
            ss >> c.center[j];
        }
        ss >> c.radialDistortion;
        // Iterate through the correpsonding SIFT file and parse all SIFT features
        // readSIFT returns a vector of SIFT features for the corresponding image/camera
        c.SiftVector = readSIFT(siftFolder + "/" + c.name.substr(0,c.name.find(".")) + ".sift",i);
        if (c.SiftVector.size() == 0){
            cerr << "Error importing SIFT features for " << c.name << "\n";
        }
        cameras[i] = c;
    }

    // Read in dense model for use in nearest neighbour calculations later
    ANNpointArray densePts = annAllocPts(maxPoints, 3);
    double eps = 0.01;
    int k = 5;  //was 5 before 
    int dim = 3;
    ANNidxArray nnIdx = new ANNidx[k];
    ANNdistArray dists = new ANNdist[k]; 
    double normals[maxPoints][dim];
    int numPoints = 0;
    for (int i = 0; i < 13; i++){
        denseStream.getline(line, 1028);
    }
    int check = 1;
    if(!denseStream.getline(line, 1028)){
        check=0;
    }
    while(numPoints < maxPoints && check) {
        istringstream ss(line);
        for (int i = 0; i < dim; i++){
            ss >> densePts[numPoints][i];
        }
        for (int i = 0; i < dim; i++){
            ss >> normals[numPoints][i];
        }
        numPoints++;
        if(!denseStream.getline(line, 1028)){
            check=0;
        }
    }
    ANNkd_tree* kdTree = new ANNkd_tree(densePts, numPoints, dim);
    denseStream.close();


    // Create sparse point with corresponding feature indices and normal plane
    // For now, normal plane is just taken directly from the nearest neighbour
    // Aggregate SIFT features for sparse point and calculate homography using
    // rotation and translation. For each camera.
    sparseStream.getline(line, 1028);
    sparseStream.getline(line, 1028);
    int numSparsePoints = stoi(line);

    vector<sparseModelPoint> sparsePoints(numSparsePoints);
    ANNpoint queryPt = annAllocPt(dim);  

    for (int i = 0; i < numSparsePoints; i++){
        sparseStream.getline(line, 1028);
        istringstream ss(line);
        sparseModelPoint smp;
        
        for (int j = 0; j < 3; j++){
            ss >> smp.point[j];
            queryPt[j] = smp.point[j];
        }
        kdTree->annkSearch(queryPt, k, nnIdx, dists, eps);  
        // for (int k2 = 0; k2 < k; k2++){
        //     for (int j = 0; j < 3; j++){
        //         smp.normal[j] += normals[nnIdx[k2]][j]; 
        //     }

        // }
        // for (int j = 0; j < 3; j++){
        //     smp.normal[j] = smp.normal[j]/k;
        // }

        for (int j = 0; j < 3; j++){
            smp.normal[j] = normals[nnIdx[0]][j];
            // smp.normal[i] = no
        } 

        // for (int j = 0; j < 3; j++){
        //     smp.normal[j] = normals[nnIdx[0]][j];
        // }

        // for (int j = 0; j < 3; j++){
        //     smp.normal[j] = 0.001;
        // }
        // smp.normal[2] = 1;

        int temp;
        for (int j = 0; j < 3; j++){
            ss >> temp;
        }

        string ns;
        ss >> ns;
        int numSift = stoi(ns);
        vector<sparseSiftFeature> sparseSifts(numSift);
        cout << "Point" << i << ": (" << smp.point[0] << "," << smp.point[1] << "," << smp.point[2] << "): \n";
        for (int j = 0; j < numSift; j++){
            sparseSiftFeature ssf;
            ss >> ns;
            int imageIndex = stoi(ns);
            ss >> ns;
            int featIndex = stoi(ns);
            Camera cam = cameras[imageIndex];
            ssf.Sift = cam.SiftVector[featIndex];
            /** Side Note: ssf.modelXY and ssf.Sfit.point both refer to 
            *      same point in the image, just with different coordinate
            *      systems. 
            *      The modelXY has the origin at the centre of the image.
            *      The pointXY obtained from SIFT has the origin in the bottom left corner
            **/
            ss >> ssf.modelXY[0];
            ss >> ssf.modelXY[1];
            // cam: list feature index and camera index for each point. It doesn't tell you which camera
            // belongs to which sift descriptor. Tian created cam to place it into camera. c is the same.
            // cam and smp are both COPIES, while &ssf is a POINTER => working on ssf edits the real cam.
            // ssf: sparse sift feature
            // smp: sparse model point
            cout << "\t" << ssf.Sift.point[0] << "," << ssf.Sift.point[1] << " in image " << cam.name << "\n"; 
            computeRotation(cam, &ssf, smp);
            computeTranslation(cam, &ssf, smp);
            computeHomography(cam, &ssf, smp.normal);
            // a P becomes a VIP
            createVIP(cam, siftFolder + "/" + cam.name, &ssf, "Point"+to_string(i)+"Sift"+to_string(j));
            sparseSifts[j] = ssf;
        }
        smp.features = sparseSifts;
        sparsePoints[i] = smp;
    }
    sparseStream.close();

    // Every point in sparsePoints should now have a list of VIPs
    // TODO: Get rid of extra information

    delete [] nnIdx;
    delete [] dists;
    delete kdTree;
    annClose();
}

// Function that obtains the coordinates of the SIFT feature as
// well as scale and orientation
vector<SiftFeature> readSIFT(string file_name, int index){
    std::string line;
    int name, version, npoint, nLocDim, nDesDim, sift_eof, sorted = 0;
    std::ifstream fs (file_name);
    if (!fs.is_open()) {
        std::cout << "Error reading file! " << strerror(errno) << "\n";
        return vector<SiftFeature>();
    }
    fs.read((char *)&name, sizeof(int));
    fs.read((char *)&version, sizeof(int));
    if(IsValidFeatureName(name) && IsValidVersionName(version)) {
        fs.read((char *)&npoint, sizeof(int));
        fs.read((char *) &nLocDim, sizeof(int));
        fs.read((char *) &nDesDim, sizeof(int));
        vector<SiftFeature> features(npoint);
        if(npoint>0 && nLocDim >0) {
            float *data = new float [nLocDim*npoint];
            char *descData = new char [nDesDim*npoint];
            
            fs.read((char *)data, nLocDim *npoint*sizeof(float));

            // Don't actually need the description here
            // fs.read((char *) descData, nDesDim*npoint*sizeof(unsigned char));

            // Here, we iterate through all SIFT features and add them to 
            // the feature vector
            for (int i = 0; i < npoint; i++){
                SiftFeature s;
                for (int j = 0; j < 3; j++){
                     s.point[j] = data[i*nLocDim + j];
                }
                s.size = data[i*nLocDim + 3];
                s.orientation = data[i*nLocDim + 4];
                // s.description = &descData[i*nDesDim];
                s.cameraIndex = index;
                features[i] = s;
            }

            fs.read((char *) &sift_eof,sizeof(int));
            fs.close();
        }else {
            fs.close();
            return vector<SiftFeature>();
        }
        // Return the feature vector.
        return features;
    }else {
        fs.close();
        return vector<SiftFeature>();
    }
}

void computeTranslation(Camera c, sparseSiftFeature *s, sparseModelPoint smp){
    // Not sure if doing this right, but essentially making the centre of the normal
    // plane a distance of focalLength away from the 3D point
    // Treating the centre of the camera as the centre of the SIFT feature. So when warping,
    // may need to crop and centre the SIFT feature first
    // Vector is from normal plane to camera plane. Easily reversed
    double c2[3];
    double distance[3];
    for (int i = 0; i < 3; i++){
        distance[i] = c.center[i] - smp.point[i];
    }
    double d = 1; // sqrt(distance[0]*distance[0]+distance[1]*distance[1] + distance[2]*distance[2]);
    for (int i = 0; i < 3; i++){
        c2[i] = smp.point[i] + d*smp.normal[i];
        // c2[i] = smp.normal[i]*smp.point[i];
    }  
    Mat C1 = Mat(3,1,CV_64FC1,c.center);
    Mat C2 = Mat(3,1,CV_64FC1,c2); 
    Mat R1 = Mat(3,3,CV_64FC1, s->R1);
    Mat t = R1*(C2-C1);

    for (int i = 0; i < 3; i++){
            s->translation[i] = t.at<double>(i,0);
        }

}

// TODO: I can't math D:
void computeRotation(Camera c, sparseSiftFeature *s, sparseModelPoint smp){
    // Turning quaternion to rotation matrix: https://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
    // Confirmation: https://groups.google.com/forum/#!topic/vsfm/V4lhITH2yHw
    double to_camera[3][3];
    double w = c.quaternion[0];
    double x = c.quaternion[1];
    double y = c.quaternion[2];
    double z = c.quaternion[3];
    to_camera[0][0] = 1 - 2*y*y - 2*z*z;
    to_camera[0][1] = 2*x*y - 2*z*w;
    to_camera[0][2] = 2*x*z + 2*y*w;
    to_camera[1][0] = 2*x*y + 2*z*w;
    to_camera[1][1] = 1 - 2*x*x - 2*z*z;
    to_camera[1][2] = 2*y*z - 2*x*w;
    to_camera[2][0] = 2*x*z - 2*y*w;
    to_camera[2][1] = 2*y*z + 2*x*w;
    to_camera[2][2] = 1 - 2*x*x - 2*y*y;

    Mat vec1;
    Mat vec2;
    Mat X = Mat(3,1,CV_64FC1,smp.normal);
    Mat up = Mat::zeros(3,1,CV_64FC1);
    up.col(0).row(2) = 1;
    vec1 = X.cross(up);
    vec1 = vec1/norm(vec1);
    vec2 = vec1.cross(X);
    vec2 = vec2/norm(vec2);

    Mat R2 = Mat::zeros(3,3,CV_64FC1);
    R2.row(0) = vec1.t();
    R2.row(1) = vec2.t();
    R2.row(2) = X.t();

    Mat R1 = Mat(3,3,CV_64FC1,to_camera);
    R1 = R1.t();
    // cout << "toCamera:" << XY_to_camera << "\n";
    // Mat trans = Mat(3,1,CV_64FC1,s->translation);
    // trans = XY_to_camera*trans;
    // Mat Normal_to_XY = MakeRotationMatrix(smp);
    // cout << "toXY:" << Normal_to_XY << "\n";
    // Mat rot = XY_to_camera*Normal_to_XY;
    // // rot = rot.inv();
    Mat R = R2*R1.t();

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++){
            s->rotation[i][j] = R.at<double>(i,j);
            s->R1[i][j] = R1.at<double>(i,j);
        }
    }

    //cout << "Rotation: " << rot << "\n";
    
}

void computeHomography(Camera c, sparseSiftFeature *s, double normal[3]){
    // Using equation 3 from the paper because that seems
    // kinda sorta right

    Mat K1 = Mat::zeros(3,3,CV_64FC1);
    K1.at<double>(0,0) = c.focalLength;
    K1.at<double>(1,1) = c.focalLength;
    K1.at<double>(2,2) = 1;

    Mat f = Mat::zeros(3,1,CV_64FC1);
    f.at<double>(0,2) = c.focalLength;

    Mat up = Mat::zeros(3,1,CV_64FC1);
    up.at<double>(0,2) = 1;
    
    Mat R = Mat(3,3,CV_64FC1, s->rotation);
    Mat R1 = Mat(3,3,CV_64FC1, s->R1);
    Mat t = Mat(3,1,CV_64FC1, s->translation);
    Mat n = Mat(3,1,CV_64FC1, normal);
    Mat c1 = Mat(3,1,CV_64FC1, c.center);

    Mat d1 = Mat::zeros(3,3,CV_64FC1);
    d1 = (R1*n).t()*R1*f;
    double d1f = d1.at<double>(0,0);

    Mat H = (R+R*t*(((R1*n).t())/d1f))*K1.inv();
    cout << "t " << t <<endl;
    cout << "H " << H <<endl;
    // Mat H = (R+R*t*n.t());

    // double tn[3][3];
    // for (int i = 0; i < 3; i++){
    //     for (int j =0; j < 3; j++){
    //         tn[i][j] = s->translation[i]*normal[j];
    //     }
    // }
    for (int i =0; i < 3; i++){
        for (int j =0; j < 3; j++){
            s->H[i][j] = H.at<double>(i,j);
        }
    }
}

void createVIP(Camera c, string imageName, sparseSiftFeature *s, string patchName){
    // Treating the size of the patch as 10*size of sift feature
    int size = s->Sift.size*10;
    Mat image, cropped;
    // Read the image
    image = imread(imageName, CV_LOAD_IMAGE_COLOR);
    // Crop the image to the SIFT patch (coordinate system from upper left corner)
    int w = image.cols;
    int h = image.rows;

    

    int x = s->Sift.point[0] - size/2;
    int y = s->Sift.point[1] - size/2;
    if (x < 0) {
        x = 0;
    }
    if (x+size >= w) {
        x = w - size - 1;
    }
    if (x < 0) {
        x = 0;
        size = w-1;
    }
    if (y < 0) {
        y = 0;
    }
    if (y+size >= h) {
        y = h - size - 1;
    }
    if (y < 0) {
        y = 0;
    }
    // Mat K1 = Mat::zeros(3,3,CV_64FC1);
    // K1.at<double>(0,0) = 1.2;
    // K1.at<double>(1,1) = 1.2;
    // K1.at<double>(2,2) = 1;
    // K1.at<double>(0,2) = x/2;
    // K1.at<double>(1,2) = y/2;
    cropped = image(Rect(x, y, size, size));

    // Warp the image
    Mat K2 = Mat::zeros(3,3,CV_64FC1);
    K2.at<double>(0,0) = c.focalLength; // c.focalLength;
    K2.at<double>(1,1) = c.focalLength; // c.focalLength;
    K2.at<double>(2,2) = 1;// homogeneous coord. 3d (world) -> 2d (pixel)
    K2.at<double>(0,2) = x/2;
    K2.at<double>(1,2) = y/2;
    Mat H = Mat(3,3,CV_64FC1,s->H);
    H = K2*H;
    // H = H.inv();
    // H = H*K1.inv();

    // Mat S = Mat::eye(3,3,CV_64F);
    // S.at<double>(0,0) = size;
    // S.at<double>(1,1) = size;
    // H = S*H*S.inv();
    // H = H.inv();
    Mat invH = H.clone().inv();
    // invH.inv();

    Mat corners = Mat(3,4,CV_64FC1);
    corners.at<double>(0,0) = 0;
    corners.at<double>(0,1) = size;
    corners.at<double>(0,2) = size;
    corners.at<double>(0,3) = 0;
    corners.at<double>(1,0) = 0;
    corners.at<double>(1,1) = 0;
    corners.at<double>(1,2) = size;
    corners.at<double>(1,3) = size;
    corners.at<double>(2,0) = 1;
    corners.at<double>(2,1) = 1;
    corners.at<double>(2,2) = 1;
    corners.at<double>(2,3) = 1;

    Mat homographyCorners = invH*corners;
    homographyCorners.row(0) = homographyCorners.row(0)/homographyCorners.row(2);
    homographyCorners.row(1) = homographyCorners.row(1)/homographyCorners.row(2);
    double minX, minY, maxX, maxY;
    minMaxLoc(homographyCorners.row(0),&minX, &maxX);
    minMaxLoc(homographyCorners.row(1),&minY, &maxY);

    int width = (int)maxX - minX;
    int height = (int)maxY - minY;
    if (width > 10000 || height > 10000){
        return;
    }
    double subX = 0;
    double subY = 0;
    if (minX < 0) {
        subX = minX;
    }
    if (minY < 0) {
        subY = minY;
    }
    Mat T = Mat::eye(3, 3, CV_64F);
    T.at<double>(0,2) = -minX;
    T.at<double>(1,2) = -minY;
    invH = T*invH;
    H = invH.clone(); // the sacred line

    int sizes[3] = {height, width, 3};
    Mat warp(3, sizes, CV_8UC(1), Scalar::all(0));
    // warpPerspective(flippedCrop, warp, H, warp.size());
    // imwrite(patchName + "VIPFlip.jpg", warp);
    warpPerspective(cropped, warp, H, warp.size());
    /// Displaying images in window
    namedWindow("Warped",CV_WINDOW_NORMAL);
    resizeWindow("Warped",400,400);
    namedWindow("Original",CV_WINDOW_NORMAL);
    moveWindow("Original",100,100);
    resizeWindow("Original",400,400);
    imshow("Original",cropped);
    imshow("Warped",warp);
    cv::waitKey(0);

    imwrite(patchName + ".jpg", cropped); // Image patches look kinda sorta right :/
    imwrite(patchName + "VIP.jpg", warp);
}

Mat MakeRotationMatrix(sparseModelPoint smp) {
    cout << "Normal (" << smp.normal[0] << "," << smp.normal[1] << "," << smp.normal[2] << "): \n";
        
    Mat R = Mat::zeros(3,3,CV_64FC1);

    Mat X = Mat(3,1,CV_64FC1,smp.normal);
    Mat Y = Mat::zeros(3,1,CV_64FC1);
    Y.col(0).row(2) = 1;
    double dotProduct = X.at<double>(2,0);
    double angle = acos(dotProduct);
    Mat cross = X.cross(Y);
    // cout << "angle: " << angle << "Cros: " << cross << "\n";
    cross = cross/norm(cross);
    // cout << "angle: " << angle << "Cros: " << cross << "\n";
    // U.col(0) = (X.col(0) + 0);
    // U.col(1) = X.cross(Y)/norm(X.cross(Y));
    // U.col(2) = (X.cross(U.col(1))+0);

    // Mat V = Mat::zeros(3,3,CV_64FC1);
    // V.col(0) = (Y.col(0) + 0);
    // V.col(1) = (U.col(1)+0);
    // V.col(2) = Y.cross(V.col(1));


    // U.t();

    double c = cos(angle);
    double s = sin(angle);
    double t = 1.0 - c;

    double x = cross.at<double>(0,0);
    double y = cross.at<double>(1,0);
    double z = cross.at<double>(2,0);

    R.at<double>(0,0) = c + x*x*t;
    R.at<double>(1,1) = c + y*y*t;
    R.at<double>(2,2) = c + z*z*t;


    double tmp1 = x*y*t;
    double tmp2 = z*s;
    R.at<double>(1,0) = tmp1 + tmp2;
    R.at<double>(0,1) = tmp1 - tmp2;
    tmp1 = x*z*t;
    tmp2 = y*s;
    R.at<double>(2,0) = tmp1 - tmp2;
    R.at<double>(0,2) = tmp1 + tmp2;    
    tmp1 = y*z*t;
    tmp2 = x*s;
    R.at<double>(2,1) = tmp1 + tmp2;
    R.at<double>(1,2) = tmp1 - tmp2;

    return R;
}







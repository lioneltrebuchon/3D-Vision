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

using namespace std;

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
      vector<SiftFeature> VIPSiftVector;
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

vector<SiftFeature> readSIFT(string filename, int index);  // ???
vector<SiftFeature> readVIPSIFT(string filename, int index, int npoint);  // ???

/******************************
* I am a file that matches SIFT descriptors
* If you don't compile me with the ann library I will be sad :( (See ann user guide)

Compile command: 
g++ matchSift.cpp -I/path/to/ann_1.1.2/include -L/path/to/ann_1.1.2/lib 


To run, something like that:
./a.out -small ../ETH-night -large ../ETH-day -sf ../ETH-night/model.nvm -df ../ETH-day/model.nvm

* Must be given 2 flags:
*       -small: the smaller photo dataset
*       -large: the larger photo dataset
*       -sf:    sparse model for the small dataset
*       -df:    sparse model for the large dataet
*******************************/
int main(int argc, char **argv)
{
    ifstream sparseStream;
    ifstream denseStream;
    string smallSiftFolder;
    string largeSiftFolder;
    int maxPoints = 1024;
    char line[1024];

    // Read in command line arguments and open the file streams
    int i = 1;
    while (i < argc) {                       
        if (!strcmp(argv[i], "-small")) {    
            smallSiftFolder = argv[++i];        
        }
        else if (!strcmp(argv[i],"-large")){
            largeSiftFolder = argv[++i];
        } 
        else if (!strcmp(argv[i], "-sf")) {   
            sparseStream.open(argv[++i], ios::in);
            if (!sparseStream) {
                cerr << "Cannot open sparse file\n";
                exit(1);
            }
        } 
        else if (!strcmp(argv[i], "-df")) {    
            denseStream.open(argv[++i], ios::in);
            if (!denseStream) {
                cerr << "Cannot open dense file\n";
                exit(1);
            }             
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
        c.SiftVector = readSIFT(smallSiftFolder + "/" + c.name.substr(0,c.name.find(".")) + ".sift",i);
        c.VIPSiftVector = readVIPSIFT(smallSiftFolder + "/" + c.name.substr(0,c.name.find(".")) + ".VIP",i,c.SiftVector.size());
        if (c.SiftVector.size() == 0){
            cerr << "Error importing SIFT features for " << c.name << "\n";
        }
        cameras[i] = c;
    }
    sparseStream.close();

    denseStream.getline(line, 1028);
    denseStream.getline(line, 1028);
    denseStream.getline(line, 1028);
    numCameras = stoi(line);
    vector<Camera> largeCameras(numCameras);
    vector<SiftFeature> SIFTS;
    vector<SiftFeature> VIPSIFTS;
    for (int i = 0; i < numCameras; i++){
        denseStream.getline(line, 1028);
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
        c.SiftVector = readSIFT(smallSiftFolder + "/" + c.name.substr(0,c.name.find(".")) + ".sift",i);
        c.VIPSiftVector = readVIPSIFT(smallSiftFolder + "/" + c.name.substr(0,c.name.find(".")) + ".VIP",i,c.SiftVector.size());
        SIFTS.insert(SIFTS.end(), c.SiftVector.begin(), c.SiftVector.end());
        VIPSIFTS.insert(VIPSIFTS.end(), c.VIPSiftVector.begin(), c.VIPSiftVector.end());
        largeCameras[i] = c;
    }

    double eps = 0.01;
    int k = 3;
    int dim = 128;
    ANNpointArray siftPts = annAllocPts(SIFTS.size(), dim);
    ANNpointArray VIPSiftPts = annAllocPts(VIPSIFTS.size(), dim);
    for (int i = 0; i < SIFTS.size(); i++){
        for (int j = 0; j < dim; j++){
            siftPts[i][j] = SIFTS[i].description[j];
        }
    }
    for (int i = 0; i < VIPSIFTS.size(); i++){
        for (int j = 0; j < dim; j++){
            VIPSiftPts[i][j] = VIPSIFTS[i].description[j];
        }
    }
    
    ANNidxArray nnIdx = new ANNidx[k];
    ANNdistArray dists = new ANNdist[k]; 
    ANNkd_tree* siftTree = new ANNkd_tree(siftPts, SIFTS.size(), dim); // ??? + where is nearest neighbor search???
    ANNkd_tree* VIPSiftTree = new ANNkd_tree(VIPSiftPts, VIPSIFTS.size(), dim); // ??? + where is nearest neighbor search???
    ANNpoint queryPt = annAllocPt(dim);  
    for (int j = 0; j < cameras.size(); j++) {
        Camera c = cameras[j];
        double aggregateDist = 0;
        double aggregateVIPDist = 0;
        for (int k = 0; k < c.SiftVector.size(); k++) {
            SiftFeature sift = c.SiftVector[k];
            for (int i = 0; i < dim; i++){
                queryPt[i] = sift.description[i];
            }
            siftTree->annkSearch(queryPt, k, nnIdx, dists, eps);
            aggregateDist += dists[0];  
        }
        for (int k = 0; k < c.VIPSiftVector.size(); k++) {
            SiftFeature sift = c.VIPSiftVector[k];
            for (int i = 0; i < dim; i++){
                queryPt[i] = sift.description[i];
            }
            VIPSiftTree->annkSearch(queryPt, k, nnIdx, dists, eps);
            aggregateVIPDist += dists[0];  
        }
        cout << "Aggregate distance for " << c.name << ": SIFT: " << aggregateDist << " VIP: " << aggregateVIPDist << "\n";
    }
    denseStream.close();

    // Every point in sparsePoints should now have a list of VIPs
    // Get rid of extra information
    delete [] nnIdx;
    delete [] dists;
    delete siftTree;
    delete VIPSiftTree;
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
        fs.read((char *) &npoint, sizeof(int));
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
                memcpy(&s.description, &descData[i*nDesDim],nDesDim);
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

// Function that obtains the coordinates of the SIFT feature as
// well as scale and orientation
vector<SiftFeature> readVIPSIFT(string file_name, int index, int npoint){
    std::string line;
    std::ifstream fs (file_name);
    if (!fs.is_open()) {
        std::cout << "Error reading file! " << strerror(errno) << "\n";
        return vector<SiftFeature>();
    }
    int j = 0;
    vector<SiftFeature> features(npoint);
    while (getline(fs, line)) {
        istringstream iss(line);
        SiftFeature s;
        iss >> s.size >> s.orientation;
        for (int i = 0; i < 128; i++){
            iss >> s.description[i];
        }
        features[j] = s;
        j++;
    }
    fs.close();
    return features;
}
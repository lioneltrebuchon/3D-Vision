#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fcntl.h>
#include <sys/stat.h> 
#include <fstream> 

enum
{
//      READ_BUFFER_SIZE = 0x100000,
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

static void printSift(float *data,char *descData, int nLocDim, int nDesDim, int nPoints){
    for (int i = 0; i < nPoints; i++){
        std::cout << "Feature " << i << ": ";
        for (int j = 0; j < nLocDim; j++){
            std::cout << data[i*nLocDim + j];
            std::cout << " ";
        }
        std::cout << "\n\n";
        if (nDesDim == 128) {
            for (int a = 0; a < 4; a++){
                for (int b = 0; b < 4; b++){
                    for (int c =0; c < 8; c++){
                        std::cout << (int)(uint32_t)(uint8_t) descData[i*nDesDim + a*4 + b*8 + c];
                        std::cout << " ";
                    }
                    std::cout << "\n";
                }
                std::cout << "\n\n";
            }
            std::cout << "--------------------------\n";
        } else {
            for (int j = 0; j < nDesDim; j++){
                std::cout << (int)(uint32_t)(uint8_t) descData[i*nDesDim + j];
                std::cout << " ";
            }
            std::cout << "\n";
        }
    }
}


// To use this function: 
// Compile using: g++ sift.cpp
// Run with: ./a.out name_of_sift_file.sift
// When run, will print out Sift features with first row containing (x,y,z,scale,orientation) 
// The following will be 4 sets of 8x4 values representing each component of the sift descriptor
int main(int argc, char* argv[]){
    if (argc != 2) {
        std::cout << "Please give a filename! \n";
        return -1;
    }
    std::string file_name = argv[1];
    std::string line;
    int name, version, npoint, nLocDim, nDesDim, sift_eof, sorted = 0;
    std::ifstream fs (file_name);
    if (!fs.is_open()) {
        std::cout << "Error reading file! " << strerror(errno) << "\n";
        return 0;
    }
    fs.read((char *)&name, sizeof(int));
    fs.read((char *)&version, sizeof(int));
    std::cout << name << " " << version;
    std::cout << "\n";

    if(IsValidFeatureName(name) && IsValidVersionName(version)) {
        //version 2 file
        fs.read((char *)&npoint, sizeof(int));
        fs.read((char *) &nLocDim, sizeof(int));
        fs.read((char *) &nDesDim, sizeof(int));
        if(npoint>0 && nLocDim >0 && nDesDim==128) {
            float *data = new float [nDesDim*npoint];
            char *descData = new char [nDesDim*npoint];
            
            fs.read((char *)data, nLocDim *npoint*sizeof(float));
            fs.read((char *) descData, nDesDim*npoint*sizeof(unsigned char));
            printSift(data, descData, nLocDim, nDesDim, npoint);

            fs.read((char *) &sift_eof,sizeof(int));
            fs.close();
        }else {
            fs.close();
            return 0;
        }
        std::cout << npoint << " " << nLocDim << " " << nDesDim;
        return 1;
    }else {
        fs.close();
        return 0;
    }
}
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fcntl.h>
#include <sys/stat.h> 
#include <fstream> 
#include "sift.h"


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
            ResizeFeatureData(npoint,nLocDim,nDesDim);
            fs.read((char *) _locData->data(), nLocDim *npoint*sizeof(float));
            fs.read((char *) _desData->data(), nDesDim*npoint*sizeof(unsigned char));
            fs.read((char *) &sift_eof,sizeof(int));
            fs.close();
            _locData->_file_version = version;
            SetUpdated();
        }else {
            ResizeFeatureData(0, 0, 0);
            fs.close();
            return 0;
        }
        std::cout << _locData;
        return 1;
    }else {
        fs.close();
        return 0;
    }
}
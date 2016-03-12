#include <iostream>
#include <iostream>
#include <fstream>
using namespace std;
#include <stdio.h>  /* defines FILENAME_MAX */
// #ifdef WINDOWS
//     #include <direct.h>
//     #define GetCurrentDir _getcwd
// #else
//     #include <unistd.h>
//     #define GetCurrentDir getcwd
//  #endif

int main () {
  ofstream myfile1;
  ofstream myfile2;
  myfile1.open ("../dense3.nvm.cmvs/00/models/option-0000.ply", std::ios_base::app);
  myfile1 << "This is the dense description.\n";
  myfile1.close();
  

  myfile2.open ("../dense3.nvm.cmvs/dense3.0.nvm", std::ios_base::app);
  myfile2 << "This is the sparse description.\n";
  myfile2.close();


  return 0;
}
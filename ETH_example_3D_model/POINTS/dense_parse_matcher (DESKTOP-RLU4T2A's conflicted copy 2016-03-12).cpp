#include <iostream>
#include <string>
#include <vector>
#include <sstream>
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
  int amount1, amount2, jump;
  myfile1.open ("../dense3.0.nvm");
  GotoLine(myfile1,3);
  myfile1 >> jump;
  GotoLine(myfile1,jump);
  myfile1 >> amount1;
  float *data1X = new float[amount1];
  float *data1Y = new float[amount1];
  float *data1Z = new float[amount1];
  std::string line;
  for (int i=0;i<amount1;i++)
  	{std::getline(myfile1,line);
  	std::istringstream stream(line);
  	stream >> data1X[i] >> data1Y[i] >> data1Y[i];
  }

myfile1.open ("../dense3.nvm.cmvs/00/models/option-0000.ply");
  while (std::getline(myfile2, line))
  	++number_of_lines;

  myfile1.close();
  myfile2.close();


  return 0;
}
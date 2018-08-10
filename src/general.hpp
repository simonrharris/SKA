#ifndef __GENERAL_HPP_INCLUDED__
#define __GENERAL_HPP__
#endif

#include <string>
#include <fstream>
#include <vector>
#include <chrono> //timing
using namespace std;

string splitFileName(const string & str);

int openFileStream(const string & fileName, ifstream & fileStream, bool verbose=true);

int fileToVector(const string & filename, vector <string> & fileargs);

void printDuration(const chrono::high_resolution_clock::time_point start);

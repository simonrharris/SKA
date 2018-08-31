#ifndef __GENERAL_HPP_INCLUDED__
#define __GENERAL_HPP__
#endif

#include <string> //std::string
#include <fstream> //std::ifstream
#include <vector> //std::vector
#include <chrono> //std::chrono
#include "gzstream.h"

extern std::string versionNumber;

std::string splitFileName(const std::string & str);

int openFileStream(const std::string & fileName, std::ifstream & fileStream, bool verbose=true);

int openGzFileStream(const std::string & fileName, igzstream & gzFileStream, bool verbose=true);

int fileToVector(const std::string & filename, std::vector < std::string > & fileargs);

void printDuration(const std::chrono::high_resolution_clock::time_point start);

int calculateMissedSNPs(float totalSiteCount, const float variantSiteCount, const int kmerSize);

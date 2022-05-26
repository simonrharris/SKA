#ifndef __SK_WEED_HPP_INCLUDED__
#define __SK_WEED_HPP__
#endif

#include <string> //std::string
#include <set> //std::set
#include <vector> //std::vector

int addKmersByFilters (std::set < std::string > & myKmerSet, const std::vector < std::string > & myFiles, const int minSamples, const int maxSamples, int & myKmerSize, const int numSamples);

int weedKmers(const std::vector < std::string > & weedfiles, const std::string & kmerfile, const float minproportion, const float maxproportion, int minsamples, int maxsamples);
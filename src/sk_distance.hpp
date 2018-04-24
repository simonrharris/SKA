#ifndef __SK_ALIGN_HPP_INCLUDED__
#define __SK_ALIGN_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector

using namespace std;



int kmerDistance(const string & distancefile, const string & clusterfile, const vector<string> & kmerfiles, const int & maxSNPS, const float & minMatched);
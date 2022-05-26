#ifndef __SK_DISTANCE_HPP_INCLUDED__
#define __SK_DISTANCE_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector


int kmerDistance(const std::string & prefix, const bool distancefile, const bool clusterfile, const std::vector < std::string > & kmerfiles, const int maxSNPS, const float minMatched, const bool includeSingletons);
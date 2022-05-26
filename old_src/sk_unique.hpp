#ifndef __SK_UNIQUE_HPP_INCLUDED__
#define __SK_UNIQUE_HPP__
#endif

#include <string> //std::string

int uniqueKmers(const std::vector < std::string > & ingroupsamples, const std::vector < std::string > & kmerfiles, const float & minproportion, const std::string & outputfile, const bool incNs);
#ifndef __SK_FIND_HPP_INCLUDED__
#define __SK_FIND_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector

int findKmersInFasta(const std::vector < std::string > & queryfiles, const std::string & reffile, const bool snponly, const bool includerepeats, const std::string & outputfile);
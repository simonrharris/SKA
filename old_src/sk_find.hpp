#ifndef __SK_FIND_HPP_INCLUDED__
#define __SK_FIND_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector
#include <unordered_map> //std::unordered_map

extern std::unordered_map < std::string, char > getGeneticCode(int translationTable);

int findKmersInFasta(const std::vector < std::string > & queryfiles, const std::string & reffile, const bool snponly, const bool includerepeats, const bool includeproduct, const std::string & outputfile);
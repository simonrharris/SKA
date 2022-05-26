#ifndef __SK_TYPE_HPP_INCLUDED__
#define __SK_TYPE_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector

void addKmerToStringMap(std::unordered_map < std::string, std::string > & myKmerMap, const std::string & myKmer, const char myBase, const std::vector < bool > & myBits, const int numSamples);

int typeKmerFile(const std::string & queryfile, const std::string & profileFile, const std::vector< std::string > & subjectfiles);
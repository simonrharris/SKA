#ifndef __SK_ALIGN_HPP_INCLUDED__
#define __SK_ALIGN_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector
#include <unordered_map> //std::unordered_map

void filterAlignment(std::unordered_map < std::string, std::string > & myKmerMap, std::vector < int > & constantBaseVector, const int numSamples, const int minrequired, const bool & variantonly, float & sitecount, float & variantsitecount);

void printAlignment(const std::string & outputfilename, const std::unordered_map < std::string, std::string > & myKmerMap, std::vector < int > & constantBaseVector, const std::vector < std::string > sampleNames, const bool variantOnly);

int alignKmers(const float & maxmissingproportion, const std::string & outputprefix, const std::vector < std::string > & kmerfiles, const bool & variantonly, const bool & printkmers, const std::vector < std::string > & sample);
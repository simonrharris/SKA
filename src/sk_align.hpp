#ifndef __SK_ALIGN_HPP_INCLUDED__
#define __SK_ALIGN_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector
#include <unordered_map> //std::unordered_map

using namespace std;

void filterAlignment(unordered_map < string, string > & myKmerMap, vector < int > & constantBaseVector, const int numSamples, const int minrequired, const bool & variantonly);

void printAlignment(const string & outputfilename, const unordered_map < string, string > & myKmerMap, vector < int > & constantBaseVector, const vector < string > sampleNames, const bool variantOnly);

int alignKmers(const float & maxmissingproportion, const string & outputfile, const vector <string> & kmerfiles, const bool & variantonly, const vector <string> & sample);
#ifndef __KMERS_HPP_INCLUDED__
#define __KMERS_HPP__
#endif

#include <unordered_map> //std::unordered_map
#include <string> //std::string
#include <array> //std::array

using namespace std;

extern char complement_table[128];

extern char base_score[128];

extern char bases[4];

void ascii_bitstring(string & mybits);

void vectorbool_from_ascii(string & myascii, vector < bool > & mybits);

int extractMiddleBase(string & kmer, char & myChar);

int lowqualitytoN(string & mysequence,const string & myquality, int & minquality);

int printKmerFile(const unordered_map < string, array < int, 8 > > & mymap, const string outputfile, const int kmersize);

int printMergedKmerFile(const unordered_map < vector < bool >, vector < string > > & mymap, const string outfileprefix, const vector < string > & mysamples, const int kmersize);

int readKmerHeader(ifstream & fileStream, int & kmersize, vector < string > & names);

int collectSampleNames(const vector < string > & files, vector < string > & names);

int getSubsample(const vector < string > & sample, const vector < string > & names, vector < bool > & include);

void addKmerToBaseArrayMap(unordered_map <string, array < int, 8 > > & kmerMap, const string & kmer, const char base, const bool firstFile);

int applyFileKmerArrayMapFilters(unordered_map <string, array < int, 8 > > & kmerMap, const int & userfilecutoff, const float & userminmaf);

int applyFinalKmerArrayMapFilters(unordered_map <string, array < int, 8 > > & kmerMap, const int & usercovcutoff, const float & userminmaf);

void reverseVectorBoolKmerMap(unordered_map < string, vector < bool > > & kmerMap, unordered_map < vector < bool >,  vector < string > > & revKmerMap);
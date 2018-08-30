#ifndef __KMERS_HPP_INCLUDED__
#define __KMERS_HPP__
#endif

#include <unordered_map> //std::unordered_map
#include <string> //std::string
#include <array> //std::array

extern char complement_table[128];

extern char base_score[128];

extern char bases[4];

void ascii_bitstring(std::string & mybits);

void vectorbool_from_ascii(std::string & myascii, std::vector < bool > & mybits);

int extractMiddleBase(std::string & kmer, char & myChar);

int lowqualitytoN(std::string & mysequence,const std::string & myquality, int & minquality);

int printKmerFile(const std::unordered_map < std::string, std::array < int, 8 > > & mymap, const std::string outputfile, const int kmersize);

int printMergedKmerFile(const std::unordered_map < std::vector < bool >, std::vector < std::string > > & mymap, const std::string outputfile, const std::vector < std::string > & mysamples, const int kmersize);

int readKmerHeader(std::ifstream & fileStream, int & kmersize, std::vector < std::string > & names);

int collectSampleNames(const std::vector < std::string > & files, std::vector < std::string > & names);

int getSubsample(const std::vector < std::string > & sample, const std::vector < std::string > & names, std::vector < bool > & include);

void addKmerToBaseArrayMap(std::unordered_map < std::string, std::array < int, 8 > > & kmerMap, const std::string & kmer, const char base, const bool firstFile);

void applyFileKmerArrayMapFilters(std::unordered_map < std::string, std::array < int, 8 > > & kmerMap, const int & userfilecutoff, const float & userminmaf);

void applyFinalKmerArrayMapFilters(std::unordered_map < std::string, std::array < int, 8 > > & kmerMap, const int & usercovcutoff, const float & userminmaf);

void reverseVectorBoolKmerMap(std::unordered_map < std::string, std::vector < bool > > & kmerMap, std::unordered_map < std::vector < bool >,  std::vector < std::string > > & revKmerMap);

void reverseStringKmerMap(std::unordered_map < std::string, std::string > & kmerMap, std::unordered_map < std::vector < bool >, std::vector < std::string > > & revKmerMap);
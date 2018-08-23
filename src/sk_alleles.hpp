#ifndef __SK_FASTA_HPP_INCLUDED__
#define __SK_FASTA_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector
#include <unordered_map> //std::unordered_map

int reverseKmerMap(std::unordered_map < std::string, std::string > & myKmerMap, std::unordered_map < std::vector < bool >, std::vector < std::string > > & myRevKmerMap, const int alleleCount);

void addKmerToStringMap(std::unordered_map < std::string, std::string > & myKmerMap, const std::string & myKmer, const char myBase, const int alleleNumber, const int alleleCount);

int allelesToKmers(const std::vector< std::string > & alleles, const long & kmerlen);

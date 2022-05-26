#ifndef __SK_REFALIGN_HPP_INCLUDED__
#define __SK_REFALIGN_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector
#include <fstream> //std::ofstream


void addKmerToBasePositionMap(std::unordered_map < std::string, std::vector <int> > & myKmerMap, const std::string & myKmer, const int myPosition);

void printVariantSites(std::ofstream & myOutfile, const int totalbases, const std::vector < std::string > & sampleNames, const std::vector < std::string > & mySequences);

int alignKmersToReference(const std::string & reference, const std::string & outputfile, const std::vector < std::string > & kmerfiles, const int & kmerlen, const bool & includeref, const bool & maprepeats, const bool & fillall, const bool & variantonly, const std::vector < std::string > & sample, const bool circularContigs);
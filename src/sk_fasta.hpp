#ifndef __SK_FASTA_HPP_INCLUDED__
#define __SK_FASTA_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector

int fastaToKmers(const std::vector < std::string > & fastas, const std::string & outfilename, const long & kmerlen, const bool circularContigs);

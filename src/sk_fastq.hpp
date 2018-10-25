#ifndef __SK_FASTQ_HPP_INCLUDED__
#define __SK_FASTQ_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector

int fastqToKmers(const std::vector < std::string > & fastqs, const std::string & outfilename, const int & kmerlen, const int & minquality, const int & userfilecutoff, const int & usercovcutoff, const float & userminmaf, const bool printAlleles);

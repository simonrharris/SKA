#ifndef __SK_FASTQ_HPP_INCLUDED__
#define __SK_FASTQ_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector

using namespace std;

int fastqToKmers(const vector<string> & fastqs, const string & outfilename, const int & kmerlen, const int & userminquality, const int & userfilecutoff, const int & usercovcutoff, const float & userminmaf);

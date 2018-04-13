#ifndef __SK_FASTQ_HPP_INCLUDED__
#define __SK_FASTQ_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector

using namespace std;

int fastqToKmers(vector<string> fastqs, string outfilename, int kmerlen, int userminquality, int userfilecutoff, int usercovcutoff, float userminmaf);

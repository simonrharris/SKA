#ifndef __SK_FASTA_HPP_INCLUDED__
#define __SK_FASTA_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector

using namespace std;

int fastaToKmers(const vector<string> & fastas, const string & outfilename, const long & kmerlen);

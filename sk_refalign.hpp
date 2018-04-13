#ifndef __SK_REFALIGN_HPP_INCLUDED__
#define __SK_REFALIGN_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector

using namespace std;


int alignKmersToReference(string reference, string outputfile, vector<string> kmerfiles, int kmerlen, bool includeref, bool maprepeats);
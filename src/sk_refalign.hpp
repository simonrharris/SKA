#ifndef __SK_REFALIGN_HPP_INCLUDED__
#define __SK_REFALIGN_HPP__
#endif

#include <string> //std::string
#include <vector> //std::vector

using namespace std;


int alignKmersToReference(const string & reference, const string & outputfile, const vector<string> & kmerfiles, const int & kmerlen, const bool & includeref, const bool & maprepeats);
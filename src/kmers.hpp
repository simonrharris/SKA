#ifndef __KMERS_HPP_INCLUDED__
#define __KMERS_HPP__
#endif

#include <unordered_map> //std::unordered_map
#include <string> //std::string
#include <array> //std::array
#include <chrono> //timing
#include <fstream> //ifstream

using namespace std;

extern char complement_table[128];

extern char base_score[128];

extern char bases[4];

char complement(const char & n);

bool reverse_is_min(const string & mystring, const int & kmerlen);

int ascii_codons(string & myDNA);

int lowqualitytoN(string & mysequence,const string & myquality, int & minquality);

int printkmerfile(unordered_map<string, array<int,5> > & mymap, string outfilename, int kmersize);

int readKmerHeader(ifstream & fileStream);
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

int fileToVector(const string & filename, vector<string> & fileargs);

char complement(const char & n);

bool reverse_is_min(const string & mystring, const int & kmerlen);

void ascii_codons(string & myDNA);

void codons_from_ascii(string & myascii);

void ascii_bitstring(string & mybits);

void vectorbool_from_ascii(string & myascii, vector < bool > & mybits);

int lowqualitytoN(string & mysequence,const string & myquality, int & minquality);

int printkmerfile(unordered_map<string, array<int,5> > & mymap, string outprefix, int kmersize);

int readKmerHeader(ifstream & fileStream, int & kmersize, vector < string > & names);

int collectSampleNames(const vector < string > & files, vector < string > & names);

int getSubsample(const vector < string > & sample, const vector < string > & names, vector < bool > & include);
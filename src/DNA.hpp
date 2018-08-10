#ifndef __DNA_HPP_INCLUDED__
#define __DNA_HPP__
#endif

#include <string> //std::string
#include <array> //std::array
#include <fstream> //ifstream
#include "gzstream.h"

using namespace std;

extern char complement_table[128];

extern char base_score[128];

extern char bases[4];

char complement(const char & n);

bool reverseComplementIsMin(const string & mystring);

bool reverseComplementIfMin(string & mystring);

int readNextFastqSequence(igzstream & gzfileStream, const string & filename, string & sequence, string & quality);

int lowqualitytoN(string & mysequence,const string & myquality, const int & userminquality, int adjustment=33);

int IUPACToN(string & mysequence);

void ascii_codons(string & myDNA);

void codons_from_ascii(string & myascii);


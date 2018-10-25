#ifndef __DNA_HPP_INCLUDED__
#define __DNA_HPP__
#endif

#include <string> //std::string
#include <array> //std::array
#include <fstream> //ifstream
#include "gzstream.h"

extern char complement_table[128];

extern char base_score[128];

extern char complement_base_score[5];

extern char bases[5];

extern char complement_bases[5];

char complement(const char & n);

bool reverseComplementIsMin(const std::string & mystring);

bool reverseComplementIfMin(std::string & mystring);

int lowqualitytoN(std::string & mysequence,const std::string & myquality, const int & userminquality, int adjustment=33);

int IUPACToN(std::string & mysequence);

int ascii_codons(std::string & myDNA);

std::string codons_from_ascii(std::string & myascii);

int readNextFastaSequence(igzstream & gzfileStream, const std::string & filename, std::string & name, std::string & sequence);

int readNextFastqSequence(igzstream & gzfileStream, const std::string & filename, std::string & sequence, std::string & quality);

int countSequencesinFasta(const std::string & filename, int & sequenceCount);

int circulariseSequence(std::string & mySequence, const int myKmerLength);


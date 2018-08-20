//g++ -O3 -std=c++0x src/kmers.cpp -lz -o bin/sk_fastq
#include <unordered_map> //std::unordered_map
#include <string> //std::string
#include <array> //std::array
#include <vector> // std::vector
#include <algorithm> //std::reverse std::transform
#include <sstream> // std::istringstream
#include <chrono> //timing
#include <cassert> //std::assert
#include <fstream> //std::ifstream
#include <iostream> //std::cout
#include <cmath> //std::pow
#include <set> //std::set
#include "gzstream.h"

using namespace std;

char complement_table[128] = {
      '-',
	  '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
	  '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
	  '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
	  '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
	  'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O','P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', 
	  '-', '-', '-', '-', '-', '-', 
	  't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o','p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 
	  '-', '-', '-', '-', '-'
};

char base_score[128] = {
      5,
	  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	  0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 
	  5, 5, 5, 5, 5, 5, 
	  0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4,
	  5, 5, 5, 5, 5
};

char bases[5] = {
          'A', 'C', 'G', 'T', 'N'
};

char complement(const char & n){   
	return complement_table[int(n)];
}

bool reverseComplementIsMin(const string & mystring){
	for (float i = 0; i < ((mystring.length())/2); ++i){
		int j = i+1;
		if (mystring[i]<complement_table[int(mystring[mystring.length()-j])]){
			return false;
		}
		else if (complement_table[int(mystring[mystring.length()-j])]<mystring[i]){
			return true;
		}
	}
	return false;
}

bool reverseComplementIfMin(string & mystring){
	if (reverseComplementIsMin(mystring)){
		reverse(mystring.begin(), mystring.end());
		transform(mystring.begin(),mystring.end(),mystring.begin(),complement);
		return true;
	}
	return false;
}

int lowqualitytoN(string & mysequence,const string & myquality, const int & userminquality, int adjustment){
	int minquality=userminquality+adjustment;
	for (string::size_type i = 0; i<mysequence.length(); ++i){
		if (myquality[i]<minquality){
			mysequence[i]='N';
		}
	}
	return 0;
}

int IUPACToN(string & mysequence){

	for (string::iterator it=mysequence.begin(); it!=mysequence.end(); ++it){
		if (base_score[int(*it)]>4){
			cout << "Unrecognised character " << *it << "in sequence" << endl;
			return 1;
		}
		else if (base_score[int(*it)]>3){
			*it='N';
		}
	}
	return 0;
}

void ascii_codons(string & myDNA){
	for (string::size_type i = 0; i<myDNA.length(); i+=3){
		assert((i+2)<myDNA.length());
		myDNA[i/3] = (base_score[int(myDNA[i])]+(base_score[int(myDNA[i+1])]*4)+(base_score[int(myDNA[i+2])]*16))+63;
	}
	myDNA.erase(myDNA.length()/3);
}

void codons_from_ascii(string & myascii){
	for (string::size_type i = 0; i<myascii.length(); ++i){
		return;
	}
}

int readNextFastqSequence(igzstream & gzfileStream, const string & filename, string & sequence, string & quality){
	string line;
	if (gzfileStream.peek()!='@'){
		cout << filename << " is not in the correct format. Expecting header line to start with @." << endl << endl;
		return 1;
	}
  	getline(gzfileStream, line);
  	getline(gzfileStream, sequence);
  	if (gzfileStream.peek()!='+'){
		cout << filename << " is not in the correct format. Expecting separator line to start with +." << endl << endl;
		return 1;
	}
  	getline(gzfileStream, line);
  	getline(gzfileStream, quality);
  	if (quality.length()!=sequence.length()){
  		cout << filename << " is not in the correct format. Sequence and quality lines must be of equal length." << endl << endl;
		return 1;
  	}

		if (not gzfileStream.good()){
			cout << filename << " is not in the correct format. Expecting a multiple of 4 lines." << endl << endl;
			return 1;
	}
	return 0;
}
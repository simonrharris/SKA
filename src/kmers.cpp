//g++ -O3 -std=c++0x src/kmers.cpp -lz -o bin/sk_fastq
#include <unordered_map> //std::unordered_map
#include <string> //std::string
#include <array> //std::array
#include <chrono> //timing
#include <cassert> //std::assert
#include <fstream> //std::ifstream
#include <iostream> //std::cout

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
          4,
	  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	  0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 
	  4, 4, 4, 4, 4, 4, 
	  0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4,
	  4, 4, 4, 4, 4
};

char bases[4] = {
          'A', 'C', 'G', 'T'
};

char complement(const char & n){   
	return complement_table[int(n)];
}

bool reverse_is_min(const string & mystring, const int & kmerlen){
	for (string::size_type i = 0; i < kmerlen; ++i){
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

int ascii_codons(string & myDNA){
	for (string::size_type i = 0; i<myDNA.length(); i+=3){
		assert((i+2)<myDNA.length());
		myDNA[i/3] = (base_score[int(myDNA[i])]+(base_score[int(myDNA[i+1])]*4)+(base_score[int(myDNA[i+2])]*16))+63;
	}
	myDNA.erase(myDNA.length()/3);
	return 0;
}

int lowqualitytoN(string & mysequence,const string & myquality, int & minquality){

	
	for (string::size_type i = 0; i<mysequence.length(); ++i){
		if (int(myquality[i])<minquality){
			mysequence[i]='N';
		}
	}
	return 0;
	
}

int printkmerfile(unordered_map<string, array<int,5> > & mymap, string outfilename, int kmersize){
	
	ofstream kmerfile(outfilename);
	kmerfile << '#' << kmersize << "\n";
	for (auto it=mymap.begin(); it!=mymap.end(); ++it){
		string kmer=it->first;
		string base;
		bool basefound=false;
		int i=0;
		for (auto it2=it->second.begin(); it2!=it->second.end(); ++it2, ++i){
		
			if (*it2>0){
				if (basefound){
					base='N';
					break;
				}
				else{
					base=bases[i];
					basefound=true;
				}
			}
		
		}
		ascii_codons(kmer);
		kmerfile << base << kmer;
		
	
	}
	kmerfile.close();
	return 0;
}


int readKmerHeader(ifstream & fileStream){
	int kmersize;
	string header;
	getline(fileStream, header);
	if (header[0]!='#'){
		cout << "Malformed kmer file\n";
		return 0;
	}
	try {
		kmersize=stoi(header.substr(1,header.length()));
	}
	catch (int e){
		cout << "An exception occurred. Exception Nr. " << e << '\n';
		return 0;
	}
	return kmersize;
}

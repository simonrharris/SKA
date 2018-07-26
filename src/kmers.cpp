//g++ -O3 -std=c++0x src/kmers.cpp -lz -o bin/sk_fastq
#include <unordered_map> //std::unordered_map
#include <string> //std::string
#include <array> //std::array
#include <vector> // std::vector
#include <sstream> // std::istringstream
#include <chrono> //timing
#include <cassert> //std::assert
#include <fstream> //std::ifstream
#include <iostream> //std::cout
#include <cmath> //std::pow
#include <set> //std::set

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

int fileToVector(const string & filename, vector<string> & fileargs){
	ifstream fileStream;
	fileStream.open(filename, ios::in);
	if (fileStream.fail()) {
		cout << "Failed to open " << filename << "\n\n";
		return 1;
	}
	string word;
	while (fileStream >> word){
		if (word.length()>500){
			cout << "Names > 500 characters are not allowed\n";
			return 0;
		}
		fileargs.push_back(word);
	}
	return 0;
}

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

void ascii_bitstring(string & mybits){
	for (string::size_type i = 0; i<mybits.length(); i+=6){
		assert((i+5)<mybits.length());
		mybits[i/6] = ((int(mybits[i])-'0')+((int(mybits[i+1])-'0')*2)+((int(mybits[i+2])-'0')*4)+((int(mybits[i+3])-'0')*8)+((int(mybits[i+4])-'0')*16)+((int(mybits[i+5])-'0')*32))+33;
	}
	mybits.erase(mybits.length()/6);
}

void vectorbool_from_ascii(string & myascii, vector < bool > & mybits){
	mybits.resize(myascii.length()*6, false);
	for (string::size_type i = 0; i<myascii.length(); ++i){
		int bits=int(myascii[i])-33;
		for (int j=5; j>=0; --j){
			int power=int(pow(2,j));
			mybits[j+(6*i)]=int(bits/power);
			bits=bits%power;
		}
	}
	return;
}

int lowqualitytoN(string & mysequence,const string & myquality, int & minquality){
	for (string::size_type i = 0; i<mysequence.length(); ++i){
		if (myquality[i]<minquality){
			mysequence[i]='N';
		}
	}
	return 0;
	
}

int printkmerfile(unordered_map<string, array<int,5> > & mymap, string outputfile, int kmersize){
	
	ofstream kmerfile(outputfile);
	kmerfile << kmersize << endl;
	kmerfile << outputfile.substr(0, outputfile.find_last_of(".")) << endl;
	kmerfile << '"';
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

int readKmerHeader(ifstream & fileStream, int & kmersize, vector < string > & names){
	
	string line;
	//getline(fileStream, line);
	//kmersize=stoi(line);
	fileStream >> kmersize;
	if (kmersize==0){
		throw 10;
	}

	string name;
	getline(fileStream, line);
	getline(fileStream, line);

	istringstream ss(line);
	while(ss >> name){
		names.push_back(name);
	}

	return 0;
}


int collectSampleNames(const vector < string > & files, vector < string > & names){
	cout << "Collecting number of samples " << endl;
	int kmersize;
	ifstream fileStream;
	for (int s = 0; s < files.size(); ++s){
		fileStream.open(files[s], ios::in);

		if (fileStream.fail()) {
			cout << "Failed to open " << files[s] << "\n" << endl;
			return 1;
		}
		int newkmersize;
		try {
			readKmerHeader(fileStream, newkmersize, names);
		}
		catch (int e){
			cout << "An exception occurred when reading file " << files[s] << ". Please check the format. Exception Nr. " << e << '\n';
			return 1;
		}
		if (s==0){
			kmersize=newkmersize;
		}

		if (newkmersize!=kmersize){
			cout << "kmer files have different kmer sizes\n" << endl;
			return 1;
		}
		fileStream.close();
	}

	int numSamples=names.size();

	cout << "Found " << numSamples << " samples in " << files.size() << " files" << endl;
	return 0;
}



int getSubsample(const vector < string > & sample, const vector < string > & names, vector < bool > & include){
	set < string > sampleSet;

	if (sample.size()>0){
		for (auto it=sample.begin(); it!=sample.end(); ++it){
			auto ret = sampleSet.insert(*it);
			if (not ret.second){
				cout << "Warning " << *it << " is in your sample file more than once" << endl;
			}
		}
		cout << sampleSet.size() << " unique sample names found in your sample file" << endl;
	}
	else {
		for (auto it=names.begin(); it!=names.end(); ++it){
			auto ret = sampleSet.insert(*it);
		}
	}

	vector < bool > included(sampleSet.size());

	

	for (auto it=names.begin(); it!=names.end(); ++it){
		auto it2 = sampleSet.find(*it);
		if (it2 != sampleSet.end()){
			if (included[distance(sampleSet.begin(), it2)]==true){
				cout << "Warning " << *it << " is in your input files more than once. Only the first occurrence will be kept" << endl;
			}
			else{
				included[distance(sampleSet.begin(), it2)]=true;
				include[distance(names.begin(), it)]=true;		
			}
		}
	}


	for (auto it=included.begin(); it!=included.end(); ++it){
		if (not *it){
			auto it2 = sampleSet.begin();
			advance(it2, distance(included.begin(), it));
			cout << "Warning " << *it2 << " is in your sample file but not in any of your input files" << endl;
		}
	}
	cout << endl;

	return 0;
}

//g++ -O3 -std=c++0x src/sk_weed.cpp -lz -o bin/sk_weed
#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <string> //std::string
#include <vector> //std::vector
#include <chrono> //timing
#include "kmers.hpp"

using namespace std;

//bin/sk_weed outputfile 

//int main(int argc, char *argv[])
int weedKmers(const vector<string> & weedfiles, const string & kmerfile, const string & outputfile)
{

	auto start = chrono::high_resolution_clock::now();
	// Create the kmer map
	unordered_map<string, string> kmerMap;
	
	cout << "Reading " << kmerfile << "\n";
	ifstream fileStream;
	fileStream.open(kmerfile, ios::in);
	if (fileStream.fail()) {
		cout << "Failed to open " << kmerfile << "\n\n";
		return 0;
	}
	int kmersize=readKmerHeader(fileStream);
	char basebuffer[1];
	char kmerbuffer[kmersize*2/3];
	while (fileStream.read(basebuffer, sizeof(basebuffer))){
		string base (basebuffer, 1);
		fileStream.read(kmerbuffer, sizeof(kmerbuffer));
		string kmer (kmerbuffer, kmersize*2/3);
		
		auto it = kmerMap.find(kmer);//check if the kmer is in the hash
		if ( it != kmerMap.end() ){//if the kmer is in the hash
			it->second=base;//make the hash a repeat
			cout << "shouldn't be here\n";
		}
		else {
			kmerMap.insert(make_pair(kmer, base));
		}
    }
	fileStream.close();
	
	cout << kmerMap.size() << " unique kmers in map\n";
	
	int weeded=0;
	int totalweedkmers=0;
	
	for (auto it = weedfiles.begin(); it != weedfiles.end(); ++it){
		cout << "Weeding kmers from " << *it << "\n";
		
		fileStream.open(*it, ios::in);
		if (fileStream.fail()) {
			cout << "Failed to open " << *it << "\n\n";
			return 0;
		}
		int weedkmersize=readKmerHeader(fileStream);
		if (weedkmersize!=kmersize){
			cout << "Files have different kmer sizes\n";
			return 0;
		}

		while (fileStream.read(basebuffer, sizeof(basebuffer))){
			string base (basebuffer, 1);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			string kmer (kmerbuffer, kmersize*2/3);
			
			auto it2 = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it2 != kmerMap.end() ){//if the kmer is in the hash
				weeded++;
				kmerMap.erase(it2);
			}
			totalweedkmers++;
				
		
    	}
		fileStream.close();
	}
	cout << totalweedkmers << " read from " << weedfiles.size() << " files\n";
	cout << "Weeded " << weeded << " kmers from " << kmerfile << "\n";
	cout << kmerMap.size() << " kmers remaining in map\n";
	
	cout << "Writing kmers to " << outputfile << "\n";
	
	ofstream weededfile(outputfile);
	weededfile << "#" << kmersize << "\n";
	
	for (auto it=kmerMap.begin(); it!=kmerMap.end(); ++it){
		
		string kmer=it->first;
		string base=it->second;
		
		weededfile << base << kmer;
		
	}
	weededfile.close();

	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done\n";
	cout << "Total time required: " << elapsed.count() << "s\n\n";
	
	return 0;
	
	
}



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
int uniqueKmers(const vector<string> & ingroupfiles, const vector<string> & outgroupfiles, const float & minproportion, const string & outputfile)
{

	auto start = chrono::high_resolution_clock::now();
	// Create the kmer map
	unordered_map<string, int> kmerMap;
	int oldkmersize=0;
	int totalingroupkmers=0;

	float maxmissing=(1.0-minproportion)*ingroupfiles.size();

	float minrequired=ingroupfiles.size()-maxmissing;

	ifstream fileStream;
	int kmersize;
	char * kmerbuffer;

	for (int s = 0; s < ingroupfiles.size(); ++s){
		cout << "Reading " << ingroupfiles[s] << "\n";
		fileStream.open(ingroupfiles[s], ios::in);

		if (fileStream.fail()) {
			cout << "Failed to open" << ingroupfiles[s] << "\n\n";
			return 0;
		}
		kmersize=readKmerHeader(fileStream);
		if (s==0){
			oldkmersize=kmersize;
		}

		char kmerbuffer[(kmersize*2/3)+1];
		while (fileStream.read(kmerbuffer, sizeof(kmerbuffer))){
			string kmer (kmerbuffer, (kmersize*2/3)+1);
			if (kmer[0]=='N'){continue;}
			auto it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				it->second++;//make the hash a repeat
			}
			else {
				if (s<=maxmissing){
					kmerMap.insert(make_pair(kmer, 1));
				}
			}
			totalingroupkmers++;
	    }
		fileStream.close();
	}

	auto it = kmerMap.begin();
	auto endIter = kmerMap.end();

	for (; it!=endIter; ){
		if (it->second<minrequired){
			kmerMap.erase(it++);
		}
		else {
			++it;
		}
	}
	
	cout << totalingroupkmers << " kmers read from " << ingroupfiles.size() << " files\n";
	cout << kmerMap.size() << " unique kmers in map\n";
	
	int weeded=0;
	int totaloutgroupkmers=0;
	
	for (auto it = outgroupfiles.begin(); it != outgroupfiles.end(); ++it){
		cout << "Excluding kmers found in " << *it << "\n";
		
		fileStream.open(*it, ios::in);
		if (fileStream.fail()) {
			cout << "Failed to open " << *it << "\n\n";
			return 0;
		}
		int filekmersize=readKmerHeader(fileStream);
		if (filekmersize!=kmersize){
			cout << "Files have different kmer sizes\n";
			return 0;
		}
		
		char kmerbuffer[(kmersize*2/3)+1];
		while (fileStream.read(kmerbuffer, sizeof(kmerbuffer))){
			string kmer (kmerbuffer, (kmersize*2/3)+1);
			
			auto it2 = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it2 != kmerMap.end() ){//if the kmer is in the hash
				weeded++;
				kmerMap.erase(it2->first);
			}
			totaloutgroupkmers++;
    	}
		fileStream.close();
	}
	cout << totaloutgroupkmers << " read from " << outgroupfiles.size() << " files\n";
	//cout << "Removed " << weeded << " kmers from " << kmerfile << "\n";
	cout << kmerMap.size() << " kmers remaining in map\n";
	
	cout << "Writing kmers to " << outputfile << "\n";
	
	ofstream outfile(outputfile);
	outfile << "#" << kmersize << "\n";
	
	for (auto it=kmerMap.begin(); it!=kmerMap.end(); ++it){
		
		string kmer=it->first;
		
		outfile << kmer;
		
	}
	outfile.close();

	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done\n";
	cout << "Total time required: " << elapsed.count() << "s\n\n";
	
	return 0;
	
	
}



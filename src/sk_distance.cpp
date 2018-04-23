//g++ -O3 -std=c++0x splitkmer.cpp -lz -o splitkmer
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "general.hpp"
#include "kmers.hpp"
#include <chrono> //timing
using namespace std;


//int main(int argc, char *argv[])
int kmerDistance(const string & outputfile, const vector<string> & kmerfiles, const bool & cluster, const int & maxSNPS, const float & minMatched)
{

	auto start = chrono::high_resolution_clock::now();
	int numfiles=kmerfiles.size();
	
	// Create the kmer map
	unordered_map<string, string> kmerMap;
	string emptySequence (numfiles , '-');
	vector< int > kmerCounts;
	int oldkmersize=0;

	ifstream fileStream;
	for (int s = 0; s < kmerfiles.size(); ++s){
		int kmercount=0;
		cout << "Reading " << kmerfiles[s] << "\n";
		fileStream.open(kmerfiles[s], ios::in);

		if (fileStream.fail()) {
			cout << "Failed to open " << kmerfiles[s] << "\n\n";
			return 0;
		}
		int kmersize=readKmerHeader(fileStream);
		if (s==0){
			oldkmersize=kmersize;
		}

		if (kmersize!=oldkmersize){
			cout << "kmer files have different kmer sizes\n\n";
			return 0;
		}

		char basebuffer[1];
		char kmerbuffer[kmersize*2/3];

		while (fileStream.read(basebuffer, sizeof(basebuffer))){
			string base (basebuffer, 1);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			string kmer (kmerbuffer, kmersize*2/3);
			kmercount++;
			
			auto it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				it->second[s]=base[0];//make the hash a repeat
			}
			else {
				string newsequence = emptySequence;
				newsequence[s] = base[0];
				kmerMap.insert(make_pair(kmer, newsequence));
			}
    	}
		fileStream.close();
		kmerCounts.push_back(kmercount);
	}
	cout << kmerMap.size() << " kmers read from " << numfiles << " files\n";
	
	cout << "Calculating distances\n";

	auto it = kmerMap.begin();
	auto endIter = kmerMap.end();
	vector< vector<int> > pairwiseSNPs;
	pairwiseSNPs.resize( numfiles , vector<int>( numfiles , 0 ) );
	vector< vector<int> > pairwiseMatches;
	pairwiseMatches.resize( numfiles , vector<int>( numfiles , 0 ) );
	vector< vector<int> > pairwiseMismatches;
	pairwiseMismatches.resize( numfiles , vector<int>( numfiles , 0 ) );
	vector< vector<int> > pairwiseNs;
	pairwiseNs.resize( numfiles , vector<int>( numfiles , 0 ) );

	for (; it!=endIter; ){
		for (int i=0; i<numfiles; ++i){
			for (int j=i+1; j<numfiles; ++j){
				if (it->second[i]=='-' && it->second[j]=='-'){
					continue;
				}
				else if (it->second[i]=='-' || it->second[j]=='-'){
					pairwiseMismatches[i][j]++;
				}
				else if (it->second[i]=='N' || it->second[j]=='N'){
					pairwiseNs[i][j]++;
				}
				else if (it->second[i]==it->second[j]){
					pairwiseMatches[i][j]++;
				}
				else {
					pairwiseSNPs[i][j]++;
				}
			}
		}
		++it;
	}

	if (cluster){
		cout << "Finding clusters\n";
		for (int i=0; i<numfiles; ++i){
			vector< int > matches;
			for (int j=i+1; j<numfiles; ++j){
				float kmercount=min(kmerCounts[i], kmerCounts[j]);
				float percentmatched = float(pairwiseMatches[i][j]+pairwiseSNPs[i][j])/kmercount;
				cout << pairwiseSNPs[i][j] << " " << percentmatched << "\n";
				if (pairwiseSNPs[i][j]<maxSNPS && percentmatched>minMatched){
					matches.push_back(j);
				}
			}
			cout << i <<": ";
			for (auto it=matches.begin(); it<matches.end(); ++it){
				cout << *it << " ";
			}
			cout <<"\n";
		}
	}
	
	cout << "Printing distances to " << outputfile << "\n";
	
	ofstream distancefile(outputfile);
	distancefile << "File 1\tFile 2\tMatches\tMismatches\tSNPs\tNs\n";
	for (int i=0; i<numfiles; ++i){
		for (int j=i+1; j<numfiles; ++j){
			distancefile << kmerfiles[i] << "\t" << kmerfiles[j] << "\t" << pairwiseMatches[i][j] << "\t" << pairwiseMismatches[i][j] << "\t" << pairwiseSNPs[i][j] << "\t" << pairwiseNs[i][j] << "\n";
		}
	}
	distancefile.close();

	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done\n";
	cout << "Total time required: " << elapsed.count() << "s\n\n";
	
	return 0;
	
	
}



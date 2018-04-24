//g++ -O3 -std=c++0x splitkmer.cpp -lz -o splitkmer
#include <unordered_map>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "general.hpp"
#include "kmers.hpp"
#include <chrono> //timing
using namespace std;


//int main(int argc, char *argv[])
int kmerDistance(const string & distancefile, const string & clusterfile, const vector<string> & kmerfiles, const int & maxSNPS, const float & minMatched)
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
	vector< vector<int> > pairwiseNs;
	pairwiseNs.resize( numfiles , vector<int>( numfiles , 0 ) );

	for (; it!=endIter; ){
		for (int i=0; i<numfiles; ++i){
			if (it->second[i]=='-'){
				continue;
			}
			for (int j=i+1; j<numfiles; ++j){
				if (it->second[j]=='-'){
					continue;
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
	kmerMap.clear();

	if (clusterfile!=""){
		ofstream clusterout(clusterfile);
		map <int, int> clusterMap;
		vector < vector <int> > clusters;
		cout << "Printing clusters to " << clusterfile << "\n";
		clusterout << "File\tCluster\n";
		for (int i=0; i<numfiles; ++i){
			vector< int > matches;
			matches.push_back(i);
			for (int j=i+1; j<numfiles; ++j){
				float kmercount=min(kmerCounts[i], kmerCounts[j]);
				float percentmatched = float(pairwiseMatches[i][j]+pairwiseSNPs[i][j])/kmercount;
				if (pairwiseSNPs[i][j]<=maxSNPS && percentmatched>minMatched){
					matches.push_back(j);
				}
			}

			int clusternum=clusters.size();
			for ( auto it=matches.begin(); it!=matches.end(); ++it){
				auto it2 = clusterMap.find(*it);//check if the match is in the hash
				if ( it2 != clusterMap.end() ){//if the match is in the hash
					if (it2->second<clusternum){
						clusternum=it2->second;
					}
				}
			}

			if (clusternum==clusters.size()){
				clusters.push_back(vector <int>());
			}

			for ( auto it=matches.begin(); it!=matches.end(); ++it){
				auto it2 = clusterMap.find(*it);//check if the match is in the hash
				if ( it2 != clusterMap.end() ){//if the match is in the hash
					if (it2->second!=clusternum){
						for ( auto it3=clusters[it2->second].begin(); it3!=clusters[it2->second].end(); ++it3){
							clusterMap[*it3]=clusternum;
							clusters[clusternum].push_back(*it3);
						}
						clusters.erase(clusters.begin()+it2->second);
					}
				}
				else{
					clusterMap.insert(make_pair(*it, clusternum));
					clusters[clusternum].push_back(*it);
				}
			}
		}
		for ( auto it=clusterMap.begin(); it!=clusterMap.end(); ++it){
			clusterout << kmerfiles[it->first] << "\t" << it->second+1 << "\n";
		}
		clusterout.close();
	}
	
	if (distancefile!=""){
		cout << "Printing distances to " << distancefile << "\n";
		
		ofstream distanceout(distancefile);
		distanceout << "File 1\tFile 2\tMatches\tMismatches\tSNPs\tNs\n";
		for (int i=0; i<numfiles; ++i){
			for (int j=i+1; j<numfiles; ++j){
				distanceout << kmerfiles[i] << "\t" << kmerfiles[j] << "\t" << pairwiseMatches[i][j] << "\t" << (kmerCounts[i]-(pairwiseMatches[i][j]+pairwiseSNPs[i][j]+pairwiseNs[i][j]))+(kmerCounts[j]-(pairwiseMatches[i][j]+pairwiseSNPs[i][j]+pairwiseNs[i][j])) << "\t" << pairwiseSNPs[i][j] << "\t" << pairwiseNs[i][j] << "\n";
			}
		}
		distanceout.close();
	}

	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done\n";
	cout << "Total time required: " << elapsed.count() << "s\n\n";
	
	return 0;
	
	
}



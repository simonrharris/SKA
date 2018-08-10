//g++ -O3 -std=c++0x splitkmer.cpp -lz -o splitkmer
#include <unordered_map>
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <bitset>
#include <cmath>       /* ceil */
#include "general.hpp"
#include "kmers.hpp"
#include "DNA.hpp"
#include <chrono> //timing
using namespace std;


//int main(int argc, char *argv[])
int kmerDistance(const string & prefix, const bool & distancefile, const bool & clusterfile, const vector<string> & kmerfiles, const int & maxSNPS, const float & minMatched)
{

	const chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();

	int numfiles=kmerfiles.size();

	vector < string > sampleNames;
	if (collectSampleNames(kmerfiles, sampleNames)!=0){return 1;}
	int numSamples=sampleNames.size();
	
	unordered_map<string, string> kmerMap;// Create the kmer map
	string emptySequence (numSamples , '-');
	int oldkmersize=0;

	vector < int > kmerCounts(numSamples, 0);

	ifstream fileStream;
	int sampleNum=0;
	for (int s = 0; s < kmerfiles.size(); ++s){
		
		if (openFileStream(kmerfiles[s], fileStream)){return 1;};

		int kmersize;
		vector < string > names;
		
		readKmerHeader(fileStream, kmersize, names);
		

		if (s==0){
			oldkmersize=kmersize;
		}

		if (kmersize!=oldkmersize){
			cout << "kmer files have different kmer sizes\n\n";
			return 0;
		}


		char basebuffer[1];
		char kmerbuffer[kmersize*2/3];
		char asciibuffer[int(ceil(float(names.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
			string asciibits (asciibuffer, sizeof(asciibuffer));
			vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools
			
			while (fileStream.peek()!='\n' && fileStream.read(basebuffer, sizeof(basebuffer))){
				string base (basebuffer, 1);
				base[0]=toupper(base[0]);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				string kmer (kmerbuffer, kmersize*2/3);
				for (int i=0; i<names.size(); ++i){ //add the base to all samples that are true in the bitset
					if (mybits[i]){
						kmerCounts[sampleNum+i]++;
					}
				}
				
				auto it = kmerMap.find(kmer);//check if the kmer is in the hash
				if ( it != kmerMap.end() ){//if the kmer is in the hash
					for (int i=0; i<names.size(); ++i){ //add the base to all samples that are true in the bitset
						if (mybits[i]){
							it->second[sampleNum+i]=base[0];//add the base to the string
						}
					}
				}
				else {
					string newsequence = emptySequence;
					for (int i=0; i<names.size(); ++i){ //add the base to all samples that are true in the bitset
						if (mybits[i]){
							newsequence[sampleNum+i] = base[0];
						}
					}
					kmerMap.insert(make_pair(kmer, newsequence));
				}
	    	}
	    	fileStream.ignore(256,'\n');//skip the end ofline character
	    }
	    sampleNum+=int(names.size()); //add the number of samples in the file to the count of total samples
		fileStream.close();
	}
	cout << kmerMap.size() << " kmers read from " << numSamples << " samples in " << numfiles << " files\n";

	vector < vector < int > > pairwiseSNPs( numSamples , vector<int>( numSamples , 0 ) );
	vector < vector < int > > pairwiseMatches( numSamples , vector<int>( numSamples , kmerMap.size() ) );
	vector < vector < int > > pairwiseNs( numSamples , vector<int>( numSamples , 0 ) );
	
	
	cout << "Calculating SNP distances" << endl;;

	vector < bitset < 1000000 > >  kmerbitvector( numSamples, 0 );
	bitset < 1000000 > mykmerbits;
	

	while (kmerMap.size()>0){
	
		auto it = kmerMap.begin();
		auto endIter = kmerMap.end();
		int j=0;
	
		for (; it!=endIter && j<1000000; ){

			vector < int > as;
			vector < int > cs;
			vector < int > gs;
			vector < int > ts;
			int ncount=0;
			
			for (int i=0; i<numSamples; ++i){
				switch (it->second[i])
				{
					case 'A':
						as.push_back(i);
						break;
					case 'C':
						cs.push_back(i);
						break;
					case 'G':
						gs.push_back(i);
						break;
					case 'T':
						ts.push_back(i);
						break;
					case '-':
					case 'N':
						kmerbitvector[i].set(j);
						ncount++;
						break;
				}
			}

			if (as.size()==numSamples || cs.size()==numSamples || gs.size()==numSamples || ts.size()==numSamples){
				kmerMap.erase(it++);
				continue;
			}
	
			int i=0;
			for (auto it=as.begin(); it!=as.end(); ++it){
				for (auto it2=cs.begin(); it2!=cs.end(); ++it2){
					if (*it<*it2){
						pairwiseSNPs[*it][*it2]++;
					}
					else {
						pairwiseSNPs[*it2][*it]++;
					}
				}
				for (auto it2=gs.begin(); it2!=gs.end(); ++it2){
					if (*it<*it2){
						pairwiseSNPs[*it][*it2]++;
					}
					else {
						pairwiseSNPs[*it2][*it]++;
					}
				}
				for (auto it2=ts.begin(); it2!=ts.end(); ++it2){
					if (*it<*it2){
						pairwiseSNPs[*it][*it2]++;
					}
					else {
						pairwiseSNPs[*it2][*it]++;
					}
				}
			}
			for (auto it=cs.begin(); it!=cs.end(); ++it){
				for (auto it2=gs.begin(); it2!=gs.end(); ++it2){
					if (*it<*it2){
						pairwiseSNPs[*it][*it2]++;
					}
					else {
						pairwiseSNPs[*it2][*it]++;
					}
				}
				for (auto it2=ts.begin(); it2!=ts.end(); ++it2){
					if (*it<*it2){
						pairwiseSNPs[*it][*it2]++;
					}
					else {
						pairwiseSNPs[*it2][*it]++;
					}
				}
			}
			for (auto it=gs.begin(); it!=gs.end(); ++it){
				for (auto it2=ts.begin(); it2!=ts.end(); ++it2){
					if (*it<*it2){
						pairwiseSNPs[*it][*it2]++;
					}
					else {
						pairwiseSNPs[*it2][*it]++;
					}
				}
			}
			if (ncount>0){
				j++;
			}
			kmerMap.erase(it++);
		}
		
		for (int i=0; i<numSamples; ++i){
			for (int j=i+1; j<numSamples; ++j){
				mykmerbits = kmerbitvector[i] | kmerbitvector[j];
				pairwiseMatches[i][j]-=mykmerbits.count();
			}
			kmerbitvector[i].reset();
		}
	}

	string dotfilename=prefix+".dot";
	ofstream dotout(dotfilename);
	dotout << "graph {" << endl;

	if (clusterfile){
		string clusterfilename=prefix+".clusters.tsv";
		ofstream clusterout(clusterfilename);
		map <int, int> clusterMap;
		vector < vector <int> > clusters;
		cout << "Printing clusters to " << clusterfilename << endl;
		clusterout << "ID\tCluster__autocolour\n";
		for (int i=0; i<numSamples; ++i){
			vector< int > matches;
			matches.push_back(i);
			//cout << i << endl;
			for (int j=i+1; j<numSamples; ++j){
				float kmercount=min(kmerCounts[i], kmerCounts[j]);
				float percentmatched = float(pairwiseMatches[i][j])/kmercount;
				if (pairwiseSNPs[i][j]<=maxSNPS && percentmatched>=minMatched){

					matches.push_back(j);
					float similarity;
					if (maxSNPS>0){
						similarity=((1.0-(float(pairwiseSNPs[i][j])/float(maxSNPS)))*3.0)+1;
					}
					else {
						similarity=1.0;
					}
					dotout << "\t" << sampleNames[i] << " -- " << sampleNames[j] << " [weight=" << similarity << "] ;" << endl;
				}
			}
			//dotout << "\t" << kmerfiles[i] << ";" << endl; //Need to find somewhere to put this where it does what I want
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
							clusters[clusternum].push_back(*it3);
						}
						clusters[it2->second].clear();
					}
				}
				else{
					clusterMap.insert(make_pair(*it, clusternum));
					clusters[clusternum].push_back(*it);
				}
			}


			for ( auto it2=clusters[clusternum].begin(); it2!=clusters[clusternum].end(); ++it2){
				auto it3 = clusterMap.find(*it2);//check if the match is in the hash
				if ( it3 != clusterMap.end() ){
					it3->second=clusternum;//segfault is on this line
				}	
			}

		}
		for ( auto it=clusterMap.begin(); it!=clusterMap.end(); ++it){
			clusterout << sampleNames[it->first] << "\t" << it->second+1 << "\n";
		}
		clusterout.close();

		int i=0;
		for ( auto it=clusters.begin(); it!=clusters.end(); ++it){
			++i;
			if (it->size()>1){
				string clusterfilename=prefix+".cluster."+to_string(i)+".txt";
				ofstream clusterout(clusterfilename);
				for ( auto it2=it->begin(); it2!=it->end(); ++it2){
					clusterout << sampleNames[*it2] << "\n";
				}
				clusterout.close();
			}
		}
	}

	dotout << "}" << endl;
	dotout.close();
	
	if (distancefile){
		string distancefilename=prefix+".distances.tsv";
		cout << "Printing distances to " << distancefilename << "\n";
		
		ofstream distanceout(distancefilename);
		distanceout << "File 1\tFile 2\tMatches\tMismatches\tSNPs\n";
		for (int i=0; i<numSamples; ++i){
			for (int j=i+1; j<numSamples; ++j){
				distanceout << sampleNames[i] << "\t" << sampleNames[j] << "\t" << pairwiseMatches[i][j] << "\t" << (kmerCounts[i]-pairwiseMatches[i][j]+kmerCounts[j]-pairwiseMatches[i][j]) << "\t" << pairwiseSNPs[i][j] << "\n";
			}
		}
		distanceout.close();
	}

	printDuration(start);
	
	return 0;
	
	
}



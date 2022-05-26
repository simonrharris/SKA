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


int kmerDistance(const std::string & prefix, const bool distancefile, const bool clusterfile, const std::vector < std::string > & kmerfiles, const int maxSNPS, const float minMatched, const bool includeSingletons)
{

	const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

	int numfiles=kmerfiles.size();

	std::vector < std::string > sampleNames;
	if (collectSampleNames(kmerfiles, sampleNames)!=0){return 1;}
	int numSamples=sampleNames.size();
	
	std::unordered_map < std::string, std::string > kmerMap;// Create the kmer map
	std::string emptySequence (numSamples , '-');
	int oldkmersize=0;

	std::vector < int > kmerCounts(numSamples, 0);

	std::ifstream fileStream;
	int sampleNum=0;
	for (int s = 0; s < kmerfiles.size(); ++s){
		
		if (openFileStream(kmerfiles[s], fileStream)){return 1;};

		int kmersize;
		std::vector < std::string > names;
		
		readKmerHeader(fileStream, kmersize, names);
		

		if (s==0){
			oldkmersize=kmersize;
		}

		if (kmersize!=oldkmersize){
			std::cout << "kmer files have different kmer sizes\n\n" << std::endl << std::endl;
			return 0;
		}


		char base;
		char kmerbuffer[kmersize*2/3];
		char asciibuffer[int(ceil(float(names.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
			std::string asciibits (asciibuffer, sizeof(asciibuffer));
			std::vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools
			
			while (fileStream.peek()!='\n' && fileStream.get(base)){
				base=toupper(base);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				std::string kmer (kmerbuffer, kmersize*2/3);
				for (int i=0; i<names.size(); ++i){ //add the base to all samples that are true in the bitset
					if (mybits[i]){
						kmerCounts[sampleNum+i]++;
					}
				}
				
				std::unordered_map < std::string, std::string >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
				if ( it != kmerMap.end() ){//if the kmer is in the hash
					for (int i=0; i<names.size(); ++i){ //add the base to all samples that are true in the bitset
						if (mybits[i]){
							it->second[sampleNum+i]=base;//add the base to the string
						}
					}
				}
				else {
					std::string newsequence = emptySequence;
					for (int i=0; i<names.size(); ++i){ //add the base to all samples that are true in the bitset
						if (mybits[i]){
							newsequence[sampleNum+i] = base;
						}
					}
					kmerMap.insert(std::make_pair(kmer, newsequence));
				}
	    	}
	    	fileStream.ignore(256,'\n');//skip the end ofline character
	    }
	    sampleNum+=int(names.size()); //add the number of samples in the file to the count of total samples
		fileStream.close();
	}
	std::cout << kmerMap.size() << " kmers read from " << numSamples << " samples in " << numfiles << " files" << std::endl;

	std::vector < std::vector < int > > pairwiseSNPs( numSamples , std::vector < int >( numSamples , 0 ) );
	std::vector < std::vector < int > > pairwiseMatches( numSamples , std::vector < int >( numSamples , kmerMap.size() ) );
	std::vector < std::vector < int > > pairwiseNs( numSamples , std::vector < int >( numSamples , 0 ) );
	
	
	std::cout << "Calculating SNP distances" << std::endl;

	std::vector < std::bitset < 1000000 > >  kmerbitvector( numSamples, 0 );
	std::bitset < 1000000 > mykmerbits;
	

	while (kmerMap.size()>0){
		
		std::vector < std::vector < int >  > baseVector (6 , std::vector < int >() ) ;
		std::unordered_map < std::string, std::string >::iterator kmit = kmerMap.begin();
		std::unordered_map < std::string, std::string >::iterator endIter = kmerMap.end();
		int j=0;
	
		for (; kmit!=endIter && j<1000000; ){

			std::fill(baseVector.begin(), baseVector.end(), ( std::vector < int >() ));

			for (int i=0; i<numSamples; ++i){

				baseVector[base_score[kmit->second[i]]].push_back(i);
				if (base_score[kmit->second[i]]==5){
					kmerbitvector[i].set(j);
				}
			}

			if (baseVector[5].size()>0){
				j++;
			}

			for (int k=0; k<4; ++k){
				for (int l=k+1; l<4; ++l){
					for (std::vector < int >::iterator it=baseVector[k].begin(); it!=baseVector[k].end(); ++it){
						for (std::vector < int >::iterator it2=baseVector[l].begin(); it2!=baseVector[l].end(); ++it2){
							//count++;
							if (*it<*it2){
								pairwiseSNPs[*it][*it2]++;
							}
							else {
								pairwiseSNPs[*it2][*it]++;
							}
						}
					}
				}
			}

			kmerMap.erase(kmit++);
		}
		
		for (int i=0; i<numSamples; ++i){
			for (int j=i+1; j<numSamples; ++j){
				mykmerbits = kmerbitvector[i] | kmerbitvector[j];
				pairwiseMatches[i][j]-=mykmerbits.count();
			}
			kmerbitvector[i].reset();
		}
	}

	std::string dotfilename=prefix+".dot";
	std::ofstream dotout(dotfilename);
	dotout << "graph {" << std::endl;

	if (clusterfile){
		std::string clusterfilename=prefix+".clusters.tsv";
		std::ofstream clusterout(clusterfilename);
		std::map <int, int> clusterMap;
		std::vector < std::vector <int> > clusters;
		std::cout << "Printing clusters to " << clusterfilename << std::endl;
		clusterout << "ID\tCluster__autocolour" << std::endl;
		for (int i=0; i<numSamples; ++i){
			std::vector < int > matches;
			matches.push_back(i);
			int matchcount=0;
			for (int j=i+1; j<numSamples; ++j){
				float kmercount=std::min(kmerCounts[i], kmerCounts[j]);
				float percentmatched = float(pairwiseMatches[i][j])/kmercount;
				if (pairwiseSNPs[i][j]<=maxSNPS && percentmatched>=minMatched){
					matchcount++;
					matches.push_back(j);
					float similarity;
					if (maxSNPS>0){
						similarity=((1.0-(float(pairwiseSNPs[i][j])/float(maxSNPS)))*3.0)+1;
					}
					else {
						similarity=1.0;
					}
					dotout << "\t" << sampleNames[i] << " -- " << sampleNames[j] << " [label=\"" << pairwiseSNPs[i][j] << "\" weight=" << similarity << "] ;" << std::endl;
				}
			}
			if (matchcount==0 && includeSingletons){
				dotout << "\t" << sampleNames[i] << " ;" << std::endl; //Need to find somewhere to put this where it does what I want
			}
			int clusternum=clusters.size();
			for ( std::vector < int >::iterator it=matches.begin(); it!=matches.end(); ++it){
				std::map <int, int >::iterator it2 = clusterMap.find(*it);//check if the match is in the hash
				if ( it2 != clusterMap.end() ){//if the match is in the hash
					if (it2->second<clusternum){
						clusternum=it2->second;
					}
				}
			}
			if (clusternum==clusters.size()){
				clusters.push_back(std::vector <int>());
			}
			for ( std::vector < int >::iterator it=matches.begin(); it!=matches.end(); ++it){
				std::map < int, int >::iterator it2 = clusterMap.find(*it);//check if the match is in the hash
				if ( it2 != clusterMap.end() ){//if the match is in the hash
					if (it2->second!=clusternum){
						for ( std::vector < int >::iterator it3=clusters[it2->second].begin(); it3!=clusters[it2->second].end(); ++it3){
							clusters[clusternum].push_back(*it3);
						}
						clusters[it2->second].clear();
					}
				}
				else{
					clusterMap.insert(std::make_pair(*it, clusternum));
					clusters[clusternum].push_back(*it);
				}
			}


			for ( std::vector < int >::iterator it2=clusters[clusternum].begin(); it2!=clusters[clusternum].end(); ++it2){
				std::map < int, int >::iterator it3 = clusterMap.find(*it2);//check if the match is in the hash
				if ( it3 != clusterMap.end() ){
					it3->second=clusternum;//segfault is on this line
				}	
			}

		}
		for ( std::map < int, int >::iterator it=clusterMap.begin(); it!=clusterMap.end(); ++it){
			clusterout << sampleNames[it->first] << "\t" << it->second+1 << "\n";
		}
		clusterout.close();

		int i=0;
		for ( std::vector < std::vector < int > >::iterator it=clusters.begin(); it!=clusters.end(); ++it){
			++i;
			if (it->size()>1){
				std::string clusterfilename=prefix+".cluster."+std::to_string(i)+".txt";
				std::ofstream clusterout(clusterfilename);
				for ( std::vector < int >::iterator it2=it->begin(); it2!=it->end(); ++it2){
					clusterout << sampleNames[*it2] << "\n";
				}
				clusterout.close();
			}
		}
	}

	dotout << "}" << std::endl;
	dotout.close();
	
	if (distancefile){
		std::string distancefilename=prefix+".distances.tsv";
		std::cout << "Printing distances to " << distancefilename << std::endl;
		
		std::ofstream distanceout(distancefilename);
		distanceout << "Sample 1\tSample 2\tMatches\tMismatches\tJaccard Index\tMash-like distance\tSNPs\tSNP distance" <<std::endl;
		for (int i=0; i<numSamples; ++i){
			for (int j=i+1; j<numSamples; ++j){
				float matches=pairwiseMatches[i][j] ;
				float mismatches=(kmerCounts[i]-pairwiseMatches[i][j]+kmerCounts[j]-pairwiseMatches[i][j]);
				float myID=matches/(mismatches+matches);
				float mashlike;
				if (myID==0){
					mashlike=1.0;
				}
				else{
					mashlike=(-1.0/((oldkmersize*2)+1))*log((2*myID)/(1+myID));
				}
				distanceout << sampleNames[i] << "\t" << sampleNames[j] << "\t" << int(matches) << "\t" << int(mismatches) << "\t" << myID << "\t" << mashlike << "\t" << pairwiseSNPs[i][j] << "\t" << pairwiseSNPs[i][j]/matches << std::endl;
			}
		}
		distanceout.close();
	}

	printDuration(start);
	
	return 0;
	
}



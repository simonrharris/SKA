#include <unordered_map>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>       // ceil
#include "kmers.hpp"
#include "general.hpp"
#include "DNA.hpp"
#include <chrono> //timing
#include "gzstream.h"


void addKmerToBasePositionMap(std::unordered_map < std::string, std::vector <int> > & myKmerMap, const std::string & myKmer, const int myPosition){
	std::unordered_map < std::string, std::vector <int> >::iterator it = myKmerMap.find(myKmer);//check if the kmer is in the hash
	if ( it != myKmerMap.end() ){//if the kmer is in the hash
		it->second.push_back(myPosition);//add the location of the match to the map
	}
	else {
		myKmerMap.insert(std::make_pair(myKmer, std::vector < int > {myPosition}));
	}
}


void printVariantSites(std::ofstream & myOutfile, const int totalbases, const std::vector < std::string > & sampleNames, const std::vector < std::string > & mySequences, const bool includeReference, const std::string & myReferenceSequence){
	std::vector < int > sitestoprint;
	sitestoprint.reserve(10000);
	std::vector < int > constantBases (4, 0);
	for (int i=0; i<totalbases; ++i){

		std::vector < int > baseVector (6, 0);

		for (int j=0; j<sampleNames.size(); ++j){
			char base=toupper(mySequences[j][i]);
			baseVector[base_score[base]]++;
		}

		int acgtCount=0;
		int constbase=-1;

		for (int i=0; i<4; ++i){
			if (baseVector[i]>0){
				acgtCount++;
				constbase=i;
			}
		}

		if (acgtCount<2){
			constantBases[constbase]++;
		}
		else {
			sitestoprint.push_back(i);
		}
		
	}

	std::cout << "Printing alignment of " << sitestoprint.size() << " variant sites" << std::endl;

	if (includeReference){
		std::string mySequence (sitestoprint.size(), '-');
		for (int j=0; j<sitestoprint.size(); ++j){
			mySequence[j]=myReferenceSequence[sitestoprint[j]];
		}
		myOutfile << mySequence << std::endl;
	}

	for (int i=0; i<sampleNames.size(); ++i){
		std::string mySequence (sitestoprint.size(), '-');
		for (int j=0; j<sitestoprint.size(); ++j){
			mySequence[j]=mySequences[i][sitestoprint[j]];
		}
		myOutfile << ">" << sampleNames[i] << std::endl << mySequence << std::endl;
	}

	std::cout << "Constant sites (a c g t):" << std::endl;
	std::cout << constantBases[0] << " " << constantBases[1]  << " " << constantBases[2]  << " " << constantBases[3]  << std::endl;
}


int alignKmersToReference(const std::string & reference, const std::string & outputfile, const std::vector < std::string > & kmerfiles, const int & kmerlen, const bool & includeref, const bool & maprepeats, const bool & fillall, const bool & variantonly, const std::vector <std::string> & sample)
{

	const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	
	// Create the kmer map
	std::unordered_map < std::string, std::vector <int> > kmerMap;
	std::set < int > revSet;
	int substringlength=(kmerlen*2)+1;
	std::string filename=splitFileName(reference);
	std::cout << "Reading " << filename << std::endl;

	int basenum=0;
	std::string refseq;
	char base;
	
	igzstream gzfileStream;

	if(openGzFileStream(reference, gzfileStream)){return 1;}//open the fasta file. May be gzipped.

	while (gzfileStream.peek() != EOF){
		std::string sequence;
		std::string name;

		if(readNextFastaSequence(gzfileStream, reference, name, sequence)){return 1;}//read the next sequence from the file

		refseq.append(sequence);
		
		int i=0;
		
		for (std::string::iterator iti = sequence.begin(), end = sequence.end()-(substringlength-1); iti != end; ++iti, ++i){
			std::string kmer=sequence.substr(i,substringlength);
			int baseposition=basenum+i+kmerlen;
			
			if (reverseComplementIfMin(kmer)){
				revSet.insert(baseposition);
			}
			extractMiddleBase(kmer, base);
			if(ascii_codons(kmer)){return 1;}
			
			addKmerToBasePositionMap(kmerMap, kmer, baseposition);

		}
		basenum+=sequence.length();
			
    }
    gzfileStream.close();


	
	std::cout << kmerMap.size() << " kmers read" << std::endl;
	
	std::vector < std::string > sampleNames;
	if (collectSampleNames(kmerfiles, sampleNames)!=0){return 1;}//get sample names from all kmerfiles

	std::vector < bool > include (sampleNames.size());
	getSubsample(sample, sampleNames, include);//get a vector of bools representing which samples to include based on a provided sample file. This also removed duplicate samples in the input files

	std::vector < std::string > includedSampleNames;
	for (std::vector < bool >::iterator it=include.begin(); it!=include.end(); ++it){//make a vector of all included sample names to help printing output files later
		if (*it){
			includedSampleNames.push_back(sampleNames[distance(include.begin(), it)]);
		}
	}

	std::cout << includedSampleNames.size() << " samples will be aligned to " << reference << std::endl;

	std::vector < std::string > sequences(includedSampleNames.size(), std::string (basenum , '-'));//create a vector to store the new sequences;

	std::ifstream fileStream;

	char kmerbuffer[kmerlen*2/3];
	int sampleNum=0;
	int includedSampleNum=0;
	std::vector < std::string > names;
	
	for (int s = 0; s < kmerfiles.size(); ++s){

		std::string filename=splitFileName(kmerfiles[s]);
		if (openFileStream(kmerfiles[s], fileStream)){return 1;};

		int kmersize;
		
		readKmerHeader(fileStream, kmersize, names);//read the header from the kmer file to get the kmer size and sample names
		
		if (kmersize!=kmerlen){
			std::cerr << std::endl << "kmer size in " << filename << " is not " << kmerlen << std::endl << std::endl;
			return 1;
		}

		std::vector < int > fileInclude;

		for (int i=0; i<names.size(); ++i){ //put the index of all sample names that are to be included into a vector
			if (include[sampleNum+i]){
				fileInclude.push_back(i);
			}
		}

		if (fileInclude.size()==0){ // if no sample names in the file are going to be included then don't read the file
			fileStream.close();
			std::cerr << "Nothing to align" << std::endl;
			sampleNum+=int(names.size());
			names.clear();
			continue;
		}

		char asciibuffer[int(ceil(float(names.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
			std::string asciibits (asciibuffer, sizeof(asciibuffer));
			std::vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools
			
			while (fileStream.peek()!='\n' && fileStream.get(base)){
				//base=toupper(base);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				std::string kmer (kmerbuffer, kmerlen*2/3);
				std::unordered_map < std::string, std::vector <int> >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
				if ( it != kmerMap.end() ){//if the kmer is in the hash
					if (it->second.size()==1 || maprepeats){//if the match is unique or the user has chosen to map repeats
						
						for (std::vector < int >::iterator itb = it->second.begin(); itb != it->second.end(); ++itb) {//iterate each match of the kmer in the reference
							char mybase=base;
							std::set < int >::iterator itc = revSet.find(*itb);//for each match see if the kmer in the reference is on the reverse strand
							if ( itc != revSet.end() ){//if so, complement the base
								mybase=complement(base);
							}

							for (int j=0; j<fileInclude.size(); ++j){ //add the base to all samples that are true in the bitset
								if (mybits[fileInclude[j]]){
									sequences[includedSampleNum+j][*itb]=mybase;//set the base of the sequence at the site matching the reference
									if (fillall || mybase!=refseq[*itb]){//if the match is a SNP or if the user has asked to map all sites in the kmer
										for (int i=*itb-kmerlen; i<=*itb+kmerlen; ++i){//for kmer sites either site of the middle base
											if (i==*itb){
												continue;//skip the middle base as we've already set that
											}
											if (sequences[includedSampleNum+j][i]=='-'){//if the base is a gap in the sample set it to the reference base in lower case
												sequences[includedSampleNum+j][i]=std::tolower(refseq[i]);
											}
										}
									}
								}
							}
						}
					}
				}
			}
			fileStream.ignore(256,'\n');//skip the end ofline character
    	}

    	for (int i=0; i<fileInclude.size(); ++i){//print the percentage of mapped bases for each sample in a file
	    	int mappedbases=basenum;
	    	int mappedNs=0;
	    	for (int j=0; j<basenum; ++j){
	    		if (sequences[includedSampleNum+i][j]=='-'){
	    			mappedbases--;
	    		}
	    	}
			std::cout << includedSampleNames[i] << ": " << float(mappedbases)/(basenum)*100 << "% of reference bases mapped" << std::endl;
    	}
    	sampleNum+=int(names.size()); //add the number of samples in the file to the count of total samples
    	includedSampleNum+=fileInclude.size();
    	names.clear();
		fileStream.close();
	}


	std::ofstream alignfile(outputfile);
	if (alignfile.fail()){
		std::cerr << std::endl << "Error: Failed to open " << outputfile << std::endl << std::endl;
		return 1;
	}

	if (includeref){
		std::cout << "Printing reference sequence to file" << std::endl;
		alignfile << ">" << filename << std::endl; 
	}
	
	if (variantonly){
		printVariantSites(alignfile, basenum, includedSampleNames, sequences, includeref, refseq);
	}
	else{
		std::cout << "Printing alignment" << std::endl;
		alignfile << refseq << std::endl;
		for (int i=0; i<includedSampleNames.size(); ++i){
			alignfile << ">" << includedSampleNames[i] << std::endl << sequences[i] << std::endl;
		}
	}

	alignfile.close();
	
	printDuration(start);

	return 0;
	
	
}



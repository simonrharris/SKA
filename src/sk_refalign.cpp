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


void printVariantSites(std::ofstream & myOutfile, const int totalbases, const std::vector < std::string > & sampleNames, const std::vector < std::string > & mySequences){
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

		if (acgtCount==1){
			constantBases[constbase]++;
		}
		else if (acgtCount>1) {
			sitestoprint.push_back(i);
		}
		
	}

	std::cout << "Printing alignment of " << sitestoprint.size() << " variant sites" << std::endl;

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


int alignKmersToReference(const std::string & reference, const std::string & outputfile, const std::vector < std::string > & kmerfiles, const int & kmerlen, const bool & includeref, const bool & maprepeats, const bool & fillall, const bool & variantonly, const std::vector <std::string> & sample, const bool circularContigs)
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
	std::vector < int > referenceContigLengths;
	
	igzstream gzfileStream;

	if(openGzFileStream(reference, gzfileStream)){return 1;}//open the fasta file. May be gzipped.

	while (gzfileStream.peek() != EOF){
		std::string sequence;
		std::string name;

		if(readNextFastaSequence(gzfileStream, reference, name, sequence)){return 1;}//read the next sequence from the file
		int seqlen=sequence.length();
		referenceContigLengths.push_back(seqlen);
		refseq.append(sequence);
		if (seqlen<substringlength){
			basenum+=seqlen;
			continue;
		}

		int baseposition=basenum;

		if (circularContigs){
			if(circulariseSequence(sequence, kmerlen)){return 1;}//if the contigs are circular, add an extra kmer length of sequence from the opposite end to each end
		}
		else {
			baseposition+=kmerlen;
		}
		
		int i=0;
		
		
		for (std::string::iterator iti = sequence.begin(), end = sequence.end()-(substringlength-1); iti != end; ++iti, ++i, ++baseposition){
			std::string kmer=sequence.substr(i,substringlength);
			
			if (reverseComplementIfMin(kmer)){
				revSet.insert(baseposition);
			}
			extractMiddleBase(kmer, base);
			if(ascii_codons(kmer)){return 1;}
			
			addKmerToBasePositionMap(kmerMap, kmer, baseposition);

		}
		basenum+=seqlen;
			
    }
    gzfileStream.close();


	
	std::cout << kmerMap.size() << " kmers read" << std::endl;
	
	std::vector < std::string > sampleNames;
	if (collectSampleNames(kmerfiles, sampleNames)!=0){return 1;}//get sample names from all kmerfiles

	std::vector < bool > include (sampleNames.size());
	getSubsample(sample, sampleNames, include);//get a vector of bools representing which samples to include based on a provided sample file. This also removed duplicate samples in the input files

	std::vector < std::string > includedSampleNames;

	if (includeref){
		includedSampleNames.push_back(filename);
	}

	for (std::vector < bool >::iterator it=include.begin(); it!=include.end(); ++it){//make a vector of all included sample names to help printing output files later
		if (*it){
			includedSampleNames.push_back(sampleNames[distance(include.begin(), it)]);
		}
	}

	std::cout << includedSampleNames.size()-int(includeref) << " samples will be aligned to " << reference << std::endl;

	std::vector < std::string > sequences(includedSampleNames.size(), std::string (basenum , '-'));//create a vector to store the new sequences;

	std::ifstream fileStream;

	char kmerbuffer[kmerlen*2/3];
	int sampleNum=0;
	int includedSampleNum=0;
	if (includeref){
		sequences[0]=refseq;
		includedSampleNum=1;
	}
	
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
								}
							}
						}
					}
				}
			}
			fileStream.ignore(256,'\n');//skip the end ofline character
    	}

    	for (int i=0; i<fileInclude.size(); ++i){//print the percentage of mapped bases for each sample in a file

    		int contigStart=0;
    		int filledto=0;

    		for (std::vector <int >::iterator rclit=referenceContigLengths.begin(); rclit!=referenceContigLengths.end(); ++rclit){//Fill bases around mapped bases or SNPs
    			filledto=contigStart;
    			for (int j=kmerlen; j<(*rclit-kmerlen); ++j){
    				if ((base_score[sequences[includedSampleNum+i][contigStart+j]]<4 && std::isupper(sequences[includedSampleNum+i][contigStart+j])) && ( fillall || sequences[includedSampleNum+i][contigStart+j] != refseq[contigStart+j])){
    					int fillfrom=std::max(std::max(contigStart+j-kmerlen, contigStart), filledto);
		    			filledto=std::min(contigStart+j+kmerlen+1, contigStart+*rclit);
		    			for (int k=fillfrom; k<filledto; ++k){
		    				if (sequences[includedSampleNum+i][k]=='-'){//if the base is a gap in the sample set it to the reference base in lower case
								sequences[includedSampleNum+i][k]=std::tolower(refseq[k]);//this needs fixing when the contigs are circular to stop it running into the next contig
							}
		    			}
    				}
    			}
    			contigStart+=*rclit;
    		}

	    	int mappedbases=basenum;
	    	int mappedNs=0;
	    	int snps=0;
	    	for (int j=0; j<basenum; ++j){
	    		if (sequences[includedSampleNum+i][j]=='-'){
	    			mappedbases--;
	    		}
	    		if (base_score[sequences[includedSampleNum+i][j]]<4 && base_score[refseq[j]]<4 && sequences[includedSampleNum+i][j]!=refseq[j] && ::isupper(sequences[includedSampleNum+i][j]) ) { 
	    			snps++;
	    		}
	    	}
			std::cout << includedSampleNames[i] << ": " << float(mappedbases)/(basenum)*100 << "% of reference bases mapped. " << snps << " SNPs found vs reference sequence. ";
			if(calculateMissedSNPs(basenum, snps, kmersize)){return 1;};
    	}
    	sampleNum+=int(names.size()); //add the number of samples in the file to the count of total samples
    	includedSampleNum+=int(fileInclude.size());
    	names.clear();
		fileStream.close();
	}


	std::ofstream alignfile(outputfile);
	if (alignfile.fail()){
		std::cerr << std::endl << "Error: Failed to open " << outputfile << std::endl << std::endl;
		return 1;
	}
	
	if (variantonly){
		printVariantSites(alignfile, basenum, includedSampleNames, sequences);
	}
	else{
		std::cout << "Printing alignment" << std::endl;
		for (int i=0; i<includedSampleNames.size(); ++i){
			alignfile << ">" << includedSampleNames[i] << std::endl << sequences[i] << std::endl;
		}
	}

	alignfile.close();
	
	printDuration(start);

	return 0;
	
	
}



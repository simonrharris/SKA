#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <string> //std::string
#include <vector> //std::vector
#include <cmath>       /* ceil */
#include "kmers.hpp"
#include "general.hpp"
#include "DNA.hpp"


int compareKmerFiles(const std::string & queryfile, const std::vector < std::string > & subjectfiles)
{

	// Create the kmer map
	std::unordered_map < std::string, char > kmerMap;

	//std::vector < std::string > sampleNames;
	//if (collectSampleNames(subjectfiles, sampleNames, false)!=0){return 1;}
	//int numSamples=sampleNames.size();
	//std::string emptySequence (numSamples , '-');
	

	int querykmersize;
	std::vector < std::string > querySampleNames;

	std::ifstream fileStream;

	if (openFileStream(queryfile, fileStream, false)){return 1;};		
		
	readKmerHeader(fileStream, querykmersize, querySampleNames);

	char base;
	char kmerbuffer[querykmersize*2/3];
	char asciibuffer[int(ceil(float(querySampleNames.size())/6))];
	while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
		std::string asciibits (asciibuffer, sizeof(asciibuffer));
		//std::vector < bool > mybits;
		//vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a verctor of bools
		while (fileStream.peek()!='\n' && fileStream.get(base)){
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			std::string kmer (kmerbuffer, querykmersize*2/3);
			
			std::unordered_map < std::string, char >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				it->second='N';
			}
			else {
				kmerMap.insert(std::make_pair(kmer, base));
			}
		}
		fileStream.ignore(256,'\n');//skip the end ofline character
	}
	querySampleNames.clear();
	fileStream.close();

	std::cout << "Subject\tKmers unique to Query\tKmers unique to Subject\tMatches\t% kmers in Query matching\t% kmers in Subject matching\tSNPs\t%ID of matching kmers\t%ID of Query kmers\t%ID of Subject kmers\tNs in Query\tNs in Subject\tNs in both" << std::endl;
	
	for (int s = 0; s < subjectfiles.size(); ++s){
		fileStream.open(subjectfiles[s], std::ios::in);
		if (fileStream.fail()) {
			std::cerr << "Failed to open " << subjectfiles[s] << std::endl << std::endl;
			return 1;
		}
		int subjectkmersize;
		std::vector < std::string > subjectSampleNames;
		
		try {
				int returnval = readKmerHeader(fileStream, subjectkmersize, subjectSampleNames);
			}
			catch (int e){
				std::cerr << "An exception occurred when reading file " << subjectfiles[s] << ". Please check the format. Exception Nr. " << e << std::endl;
				return 1;
			}
		if (subjectkmersize!=querykmersize){
			std::cout << subjectfiles[s] << " has different kmer size to query" << std::endl;
			return 1;
		}

		
		std::vector < int > kmerjustinsubject (subjectSampleNames.size(), 0);
		std::vector < int > snps (subjectSampleNames.size(), 0);
		std::vector < int > ninsubject (subjectSampleNames.size(), 0);
		std::vector < int > ninquery (subjectSampleNames.size(), 0);
		std::vector < int > ninboth (subjectSampleNames.size(), 0);
		std::vector < int > matches (subjectSampleNames.size(), 0);

		char base;
		char kmerbuffer[querykmersize*2/3];
		char asciibuffer[int(ceil(float(subjectSampleNames.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			std::string asciibits (asciibuffer, sizeof(asciibuffer));
			std::vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a verctor of bools

			int nb=0;
			int nq=0;
			int ns=0;
			int m=0;
			int sn=0;
			int js=0;

			while (fileStream.peek()!='\n' && fileStream.get(base)){
				base=toupper(base);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				std::string kmer (kmerbuffer, querykmersize*2/3);

				std::unordered_map < std::string, char >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
				if ( it != kmerMap.end() ){//if the kmer is in the hash
					if (it->second=='N' && base=='N'){
						nb++;
					}
					else if (it->second=='N'){
						nq++;
					}
					else if (base=='N'){
						ns++;
					}
					else if (it->second==base){
						m++;
					}
					else {
						sn++;
					}
				}
				else {
					js++;
				}
	    	}
	    	for (int i=0; i<subjectSampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
				if (mybits[i]){
					ninboth[i]+=nb;
					ninquery[i]+=nq;
					ninsubject[i]+=ns;
					matches[i]+=m;
					snps[i]+=sn;
					kmerjustinsubject[i]+=js;
				}
			}
	    	fileStream.ignore(256,'\n');//skip the end ofline character
	    }
		fileStream.close();
		
		for (int i=0; i<subjectSampleNames.size(); ++i){

			int kmerjustinquery=kmerMap.size()-(matches[i]+snps[i]+ninquery[i]+ninboth[i]+ninsubject[i]);
			float percentmatchquery=float(matches[i]+snps[i]+ninquery[i]+ninboth[i])/(kmerjustinquery+matches[i]+snps[i]+ninquery[i]+ninboth[i])*100;
			float percentmatchsubject=float(matches[i]+snps[i]+ninsubject[i]+ninboth[i])/(kmerjustinsubject[i]+matches[i]+snps[i]+ninsubject[i]+ninboth[i])*100;
			float percentidofmatches=float(matches[i])/(matches[i]+snps[i])*100;
			float percentidofquery=(percentidofmatches*percentmatchquery)/100;
			float percentidofsubject=(percentidofmatches*percentmatchsubject)/100;
			
			std::cout << subjectSampleNames[i] << "\t" << kmerjustinquery << "\t" << kmerjustinsubject[i] << "\t" << (matches[i]+snps[i]+ninquery[i]+ninboth[i]+ninsubject[i]) << "\t" << percentmatchquery << "\t" << percentmatchsubject << "\t" << snps[i] << "\t" << percentidofmatches << "\t" << percentidofquery << "\t" << percentidofsubject << "\t" << ninquery[i] << "\t" << ninsubject[i] << "\t" << ninboth[i] << std::endl;

		}
	}
	return 0;
	
}



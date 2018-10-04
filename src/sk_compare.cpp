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
	std::unordered_map < std::string, std::string > kmerMap;

	std::vector < std::string > sampleNames;
	if (collectSampleNames(subjectfiles, sampleNames, false)!=0){return 1;}
	int numSamples=sampleNames.size();
	std::string emptySequence (numSamples , '-');
	std::vector < int > queryKmers (numSamples, 0);

	int sampleNum=0;
	int subjectkmersize;
	int oldsubjectkmersize;
	std::vector < std::string > subjectSampleNames;

	std::ifstream fileStream;

	for (int s = 0; s < subjectfiles.size(); ++s){
		
		if (openFileStream(subjectfiles[s], fileStream, false)){return 1;};		
		
		readKmerHeader(fileStream, subjectkmersize, subjectSampleNames);

		if (s==0){
			oldsubjectkmersize=subjectkmersize;
		}

		if (subjectkmersize!=oldsubjectkmersize){
			std::cerr << "kmer files have different kmer sizes" << std::endl <<std::endl;
			return 1;
		}
		char base;
		char kmerbuffer[subjectkmersize*2/3];
		char asciibuffer[int(ceil(float(subjectSampleNames.size())/6))];
		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			std::string asciibits (asciibuffer, sizeof(asciibuffer));
			std::vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a verctor of bools
			while (fileStream.peek()!='\n' && fileStream.get(base)){
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				std::string kmer (kmerbuffer, subjectkmersize*2/3);
				
				std::unordered_map < std::string, std::string >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
				if ( it != kmerMap.end() ){//if the kmer is in the hash
					for (int i=0; i<subjectSampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
						if (mybits[i]){
							it->second[i+sampleNum]=base;
							queryKmers[i+sampleNum]++;
						}
					}
				}
				else {
					std::string newsequence = emptySequence;
					for (int i=0; i<subjectSampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
						if (mybits[i]){
							newsequence[i+sampleNum]=base;
							queryKmers[i+sampleNum]++;
						}
					}
					kmerMap.insert(std::make_pair(kmer, newsequence));
				}
			}
			fileStream.ignore(256,'\n');//skip the end ofline character
		}
		sampleNum+=int(subjectSampleNames.size()); //add the number of samples in the file to the count of total samples
		subjectSampleNames.clear();
		fileStream.close();
	}
	
	fileStream.open(queryfile, std::ios::in);
	if (fileStream.fail()) {
		std::cerr << "Failed to open " << queryfile << std::endl << std::endl;
		return 1;
	}
	int querykmersize;
	std::vector < std::string > querySampleNames;
	if (querySampleNames.size()>1){
		std::cerr << "Error: Query can only be a single sample, not a kmerge" << std::endl;
		return 1;
	}
	
	try {
			int returnval = readKmerHeader(fileStream, querykmersize, querySampleNames);
		}
		catch (int e){
			std::cerr << "An exception occurred when reading file " << queryfile << ". Please check the format. Exception Nr. " << e << std::endl;
			return 1;
		}
	if (subjectkmersize!=querykmersize){
		std::cout << queryfile << " has different kmer size to query" << std::endl;
		return 1;
	}

	std::cout << "Subject\tKmers unique to Query\tKmers unique to Subject\tMatches\t% kmers in Query matching\t% kmers in Subject matching\tSNPs\t%ID of matching kmers\t%ID of Query kmers\t%ID of Subject kmers\tNs in Query\tNs in Subject\tNs in both" << std::endl;
	
	std::vector < int > kmerjustinb (sampleNames.size(), 0);
	std::vector < int > snps (sampleNames.size(), 0);
	std::vector < int > nina (sampleNames.size(), 0);
	std::vector < int > ninb (sampleNames.size(), 0);
	std::vector < int > ninboth (sampleNames.size(), 0);
	std::vector < int > matches (sampleNames.size(), 0);

	char base;
	char kmerbuffer[querykmersize*2/3];
	char asciibuffer[int(ceil(float(querySampleNames.size())/6))];

	while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
		std::string asciibits (asciibuffer, sizeof(asciibuffer));
		while (fileStream.peek()!='\n' && fileStream.get(base)){
			base=toupper(base);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			std::string kmer (kmerbuffer, querykmersize*2/3);

			std::unordered_map < std::string, std::string >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				for (int i=0; i<sampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
					if (it->second[i]=='N' && base=='N'){
						ninboth[i]++;
					}
					else if (it->second[i]=='-'){
						kmerjustinb[i]++;
					}
					else if (it->second[i]=='N'){
						nina[i]++;
					}
					else if (base=='N'){
						ninb[i]++;
					}
					else if (it->second[i]==base){
						matches[i]++;
					}
					else {
						snps[i]++;
					}
				}
			}
			else {
				for (int i=0; i<sampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
					kmerjustinb[i]++;
				}
			}
    	}
    	fileStream.ignore(256,'\n');//skip the end ofline character
    }
	fileStream.close();
	
	for (int i=0; i<sampleNames.size(); ++i){

		int kmerjustina=queryKmers[i]-(matches[i]+snps[i]+nina[i]+ninboth[i]+ninb[i]);
		float percentmatcha=float(matches[i]+snps[i]+nina[i]+ninboth[i])/(kmerjustina+matches[i]+snps[i]+nina[i]+ninboth[i])*100;
		float percentmatchb=float(matches[i]+snps[i]+ninb[i]+ninboth[i])/(kmerjustinb[i]+matches[i]+snps[i]+ninb[i]+ninboth[i])*100;
		float percentidofmatches=float(matches[i])/(matches[i]+snps[i])*100;
		float percentidofquery=(percentidofmatches*percentmatcha)/100;
		float percentidofsubject=(percentidofmatches*percentmatchb)/100;
		
		std::cout << sampleNames[i] << "\t" << kmerjustina << "\t" << kmerjustinb[i] << "\t" << matches[i] << "\t" << percentmatcha << "\t" << percentmatchb << "\t" << snps[i] << "\t" << percentidofmatches << "\t" << percentidofquery << "\t" << percentidofsubject << "\t" << nina[i] << "\t" << ninb[i] << "\t" << ninboth[i] << std::endl;

	}
	return 0;
	
}



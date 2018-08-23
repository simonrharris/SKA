#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <string> //std::string
#include <vector> //std::vector
#include <chrono> //std::chrono
#include <cmath> //std::ceil
#include "kmers.hpp"
#include "general.hpp"
#include "DNA.hpp"

int weedKmers(const std::vector < std::string > & weedfiles, const std::string & kmerfile)
{

	const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	// Create the kmer map
	std::unordered_map < std::string, char > kmerMap;
	
	std::ifstream fileStream;

	if (openFileStream(kmerfile, fileStream)){return 1;};

	int kmersize;
	std::vector < std::string > weednames;
	
	readKmerHeader(fileStream, kmersize, weednames);
	
	char base;
	char kmerbuffer[kmersize*2/3];
	char asciibuffer[int(ceil(float(weednames.size())/6))];
	while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
		std::string asciibits (asciibuffer, sizeof(asciibuffer));
		
		while (fileStream.peek()!='\n' && fileStream.get(base)){
			base=toupper(base);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			std::string kmer (kmerbuffer, sizeof(kmerbuffer));
			kmerMap.insert(std::make_pair(kmer, base));
		}
		fileStream.ignore(256,'\n');//skip the end ofline character
    }
	fileStream.close();
	
	std::cout << kmerMap.size() << " unique kmers in map" << std::endl;
	
	int weeded=0;
	int totalweedkmers=0;
	
	for (std::vector < std::string >::const_iterator it = weedfiles.begin(); it != weedfiles.end(); ++it){

		if (openFileStream(*it, fileStream)){return 1;};

		int weedkmersize;
		std::vector < std::string > names;
		
		readKmerHeader(fileStream, weedkmersize, names);
		
		if (weedkmersize!=kmersize){
			std::cerr << "Files have different kmer sizes" << std::endl << std::endl;
			return 1;
		}
		std::string filename=std::string(*it);
		std::string outputfile=filename.substr(0, filename.find_last_of("."))+".weeded"+filename.substr(filename.find_last_of("."), filename.length());

		std::cout << "Writing weeded kmers to " << outputfile << std::endl;
	
		std::ofstream weededfile(outputfile);
		weededfile << kmersize << std::endl;
		for ( std::vector < std::string >::iterator it=names.begin(); it!=names.end(); ++it){ //print each sample name to output file stream
			weededfile << *it << " "; 
		}
		weededfile << std::endl;

		bool firstkmer=true;

		char asciibuffer[int(ceil(float(names.size())/6))];
		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
			std::string asciibits (asciibuffer, sizeof(asciibuffer));
			
			while (fileStream.peek()!='\n' && fileStream.get(base)){
				base=toupper(base);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				std::string kmer (kmerbuffer, kmersize*2/3);
				
				std::unordered_map < std::string, char >::iterator it2 = kmerMap.find(kmer);//check if the kmer is in the hash
				if ( it2 == kmerMap.end() ){//if the kmer is not in the hash
					if (firstkmer){//if it's the first kmer, print the ascii bitstring
						weededfile << asciibits;
						firstkmer=false;
					}
					weededfile << base << kmer;//print the kmer
					weeded++;
				}
				totalweedkmers++;
			}
			if (firstkmer==false){
				weededfile << std::endl;
			}
			firstkmer=true;
			fileStream.ignore(256,'\n');//skip the end ofline character
    	}
		fileStream.close();
		weededfile.close();
	}
	std::cout << totalweedkmers << " kmers prior to weeding" << std::endl;
	std::cout << weeded << " kmers remaining after weeding" << std::endl;
	std::cout << totalweedkmers-weeded << " kmers weeded (Note this may be more than the number of kmers in the weed file due to mulitple variants of a kmer being in the file)" << std::endl;

	printDuration(start);
	
	return 0;
	
	
}



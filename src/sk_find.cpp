#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <sstream> //std::stringstream
#include <string> //std::string
#include <set> //std::set
#include <vector> //std::vector
#include <cmath> //std::ceil
#include <algorithm> //std::count
#include <chrono> //std::chrono
#include <ctime> //std::localtime
#include "kmers.hpp"
#include "general.hpp"
#include "DNA.hpp"
#include "gzstream.h"


int findKmersInFasta(const std::vector < std::string > & queryfiles, const std::string & reffile, const bool snponly, const bool includerepeats, const std::string & outputfile)
{

	const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();//start the clock

	int kmersize;
	std:: vector < std::string > sampleNames;

	std::ifstream fileStream;

	if (openFileStream(queryfiles[0], fileStream, false)){return 1;};	
		
	readKmerHeader(fileStream, kmersize, sampleNames);

	std::set < std::string > fastaKmerSet;
	std::set < std::string > repeatKmerSet;

	//std::cout << "Reading " << reffile << std::endl;
	igzstream gzfileStream;

	if(openGzFileStream(reffile, gzfileStream)){return 1;}//open the fasta file. May be gzipped.

	int substringlength=(kmersize*2)+1;
	while (gzfileStream.peek() != EOF){

		while (gzfileStream.peek()=='\n'){
			gzfileStream.ignore(256,'\n');//skip the end ofline character
		}

		if (gzfileStream.peek()=='#'){
			std::string gffline;
			while (std::getline(gzfileStream,gffline)){
				if (gffline.length()>1 && gffline[0]=='#' && gffline[1]=='#'){
					std::cout << "Comment: " << gffline << std::endl;
					if (gffline=="##FASTA"){
						break;
					}
				}
			}
		}
		else if (gzfileStream.peek()!='>'){
			std::cerr << "Error: " << reffile << " is not in the correct format. Expecting header line to start with # or >." << std::endl << std::endl;
			return 1;
		}

		if (gzfileStream.peek()==EOF){
			break;
		}

		std::string sequence;
		std::string name;

		if(readNextFastaSequence(gzfileStream, reffile, name, sequence)){return 1;}//read the next sequence from the file

		if(IUPACToN(sequence)){return 1;}//convert any IUPAC characters to N and reject any unrecognised characters
		
		if (sequence.length() < substringlength){//if the subsequence is too short continue
			continue;
		}
		
		int i;

		for (i=0; i < sequence.length() - (substringlength-1); ++i){//for each base in the subsequence

			std::string kmer = sequence.substr(i,substringlength);//extract the kmer
			
			reverseComplementIfMin(kmer);//if the revcomp of the kmer is 'smaller' then revcomp it
			char base;
			if(extractMiddleBase(kmer, base)){return 1;}//extract the middle base
			if(ascii_codons(kmer)){return 1;}

			//std::cout << kmerMap[kmer][0] << std::endl;

			std::set < std::string >::iterator it = fastaKmerSet.find(kmer);//check if the kmer is in the hash
			if ( it != fastaKmerSet.end() ){//if the kmer is in the hash
				repeatKmerSet.insert(kmer);
			}
			else {
				fastaKmerSet.insert(kmer);
			}
		}
	}
	gzfileStream.close();//close the file

	// Create the kmer map
	std::unordered_map < std::string, std::array < int, 8 > > kmerMap;

	int sampleNum=0;

	int oldkmersize;
	bool firstFile=true;
	

	

	for (std::vector < std::string >::const_iterator queryit=queryfiles.begin(); queryit!=queryfiles.end(); ++queryit){//for each fasta file
		
		std::string queryfile = *queryit;
		if (openFileStream(queryfile, fileStream, false)){return 1;};		
		
		std:: vector < std::string > sampleNames;
		
		readKmerHeader(fileStream, kmersize, sampleNames);

		if (firstFile){
			oldkmersize=kmersize;
		}
		else if (kmersize!=oldkmersize){
			std::cerr << "kmer files have different kmer sizes" << std::endl <<std::endl;
			return 1;
		}

		char base;
		char kmerbuffer[kmersize*2/3];
		char asciibuffer[int(ceil(float(sampleNames.size())/6))];
		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			std::string asciibits (asciibuffer, sizeof(asciibuffer));
			std::vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a verctor of bools
			while (fileStream.peek()!='\n' && fileStream.get(base)){
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				std::string kmer (kmerbuffer, kmersize*2/3);
				
				std::unordered_map < std::string, std::array < int, 8 > >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
				if ( it != kmerMap.end() ){//if the kmer is in the hash
					it->second[base_score[base]]+=std::count(mybits.begin(), mybits.end(), true);//increment the count of the base for the kmer
				}
				else {//if the kmers isn't in the hash and we are adding from the first fastq file
					std::pair < std::unordered_map < std::string, std::array < int, 8 > >::iterator,bool>  ret = kmerMap.insert(std::make_pair(kmer, std::array < int, 8 > ()));//insert the kmer into the hash
					ret.first->second[base_score[base]]=std::count(mybits.begin(), mybits.end(), true);//increment the count of the base for the kmer
				}
			}
			fileStream.ignore(256,'\n');//skip the end ofline character
		}
		firstFile=false;
		fileStream.close();
	}





	//std::cout << repeatKmerSet.size() << " repeat kmers in set" << std::endl;
	
	//std::cout << kmerMap.size() << " unique kmers in map" << std::endl;
	
	std::unordered_map < std::string, std::vector < int > > matchMap;

	//std::cout << "Reading " << reffile << std::endl;
	//igzstream gzfileStream;

	std::ofstream vcffile(outputfile);
	vcffile << "##fileformat=VCFv4.3" << std::endl;
	std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    time_t in_time_t = std::chrono::system_clock::to_time_t(now);
    char timestr[100];
    if (std::strftime(timestr, sizeof(timestr), "%Y%m%d", std::localtime(&in_time_t))) {
        vcffile << "##fileDate=" << timestr << std::endl;
    }
	vcffile << "##source=SKA v"<< versionNumber << std::endl;
	vcffile << "##reference="<< reffile << std::endl;
	vcffile << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Matching Kmer\">" << std::endl;
	if (includerepeats){
		vcffile << "##INFO=<ID=RR,Number=0,Type=Flag,Description=\"Repeat Region\">" << std::endl;
	}
	vcffile << "##FORMAT=<ID=NS4,Number=4,Type=Integer,Description=\"Number Samples With A, C, G or T\">" << std::endl;
	vcffile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;

	int totalbases=0;

	if(openGzFileStream(reffile, gzfileStream)){return 1;}//open the fasta file. May be gzipped.

	std::unordered_map < std::string, std::vector < int > > contigAnnotation;
	std::vector < std::vector < std::string > > features;

	while (gzfileStream.peek() != EOF){

		while (gzfileStream.peek()=='\n'){
			gzfileStream.ignore(256,'\n');//skip the end ofline character
		}

		if (gzfileStream.peek()=='#'){
			std::cout << "Parsing Annotation" << std::endl;
			std::string gffline;
			while (std::getline(gzfileStream,gffline)){
				if ((gffline.length()==0) ||  (gffline[0]=='#' && (gffline.length()==1 || gffline[0]!='#'))){
					continue;
				}
				if (gffline.length()>1 && gffline[0]=='#' && gffline[1]=='#'){
					//std::cout << "Comment: " << gffline << std::endl;
					if (gffline=="##FASTA"){
						break;
					}
				}
				else {
					std::stringstream linestream;
					linestream << gffline;//convert the sequence to stringstream
					std::string word;
					std::vector < std::string > words;
					words.reserve(7);
					while (std::getline(linestream, word, '\t')){
						words.push_back(std::string(word));
					}
					std::unordered_map < std::string, std::vector < int > >::iterator it = contigAnnotation.find(words[0]);//check if the contig is in the hash
					if ( it != contigAnnotation.end() ){//if the contig is in the hash
						if (std::stoi(words[4])>it->second.size()){
							it->second.reserve(std::stoi(words[4]));
						}
						for (int i=std::stoi(words[3])-1; i<std::stoi(words[4]); ++i){
							it->second[i]=features.size();
						}
						std::cout << it->second.size() << std::endl;
					}
					else {//if the kmers isn't in the hash and we are adding from the first fastq file
						std::pair < std::unordered_map < std::string, std::vector < int> >::iterator,bool>  ret = contigAnnotation.insert(std::make_pair(words[0], std::vector < int> (10000000, -1)));//insert the kmer into the hash
						
						for (int i=std::stoi(words[3])-1; i<std::stoi(words[4]); ++i){
							ret.first->second[i]=features.size();
						}
					}
					features.push_back(words);
				}
			}
		}
		else if (gzfileStream.peek()!='>'){
			std::cerr << "Error: " << reffile << " is not in the correct format. Expecting header line to start with # or >." << std::endl << std::endl;
			return 1;
		}

		if (gzfileStream.peek()==EOF){
			break;
		}

		std::string sequence;
		std::string name;
		std::cout << features.size() << " features found in " << contigAnnotation.size() << " contigs" << std::endl;

		for (std::unordered_map < std::string, std::vector < int > >::iterator cait=contigAnnotation.begin(); cait!=contigAnnotation.end(); ++cait){
			std::cout << cait->first;
			for (std::vector < int >::iterator viit=cait->second.begin(); viit!=cait->second.end(); ++viit){
				std::cout << " " << *viit << " " << features[*viit][8] << std::endl;
			}
		}

		std::cout << "Parsing Sequence" << std::endl;
		if(readNextFastaSequence(gzfileStream, reffile, name, sequence)){return 1;}//read the next sequence from the file

		if(IUPACToN(sequence)){return 1;}//convert any IUPAC characters to N and reject any unrecognised characters

		if (sequence.length() < substringlength){//if the subsequence is too short continue
			continue;
		}

		std::unordered_map < std::string, std::vector < int > >::iterator cait = contigAnnotation.find(name);//check if the contig is in the hash
		

		for (int i=0; i < sequence.length() - (substringlength-1); ++i){//for each base in the subsequence

			std::string kmer = sequence.substr(i,substringlength);//extract the kmer
			
			bool isrev=reverseComplementIfMin(kmer);//if the revcomp of the kmer is 'smaller' then revcomp it
			char base;
			if(extractMiddleBase(kmer, base)){return 1;}//extract the middle base
			if(ascii_codons(kmer)){return 1;}

			std::unordered_map < std::string, std::array < int, 8 > >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash

			if ( it != kmerMap.end() ){//if the kmer is in the hash
				int basecount=0;
				int samplesmatching=0;
				int baseintrep=90000;
				std::string altstring="";
				std::string NS4string=";NS4=";
				//it->second[base_score[base]]++;
				int j;
				if (isrev){
					j=3;
				}
				else{
					j=0;
				}
				for (; j<4 && j>-1;){
					if (it->second[j]>0 || (j==base_score[base])){
						basecount++;
						samplesmatching+=it->second[j];
						if (bases[j]!=base){
							if (altstring.length()>0){
								altstring+=",";
							}
							if (isrev){
								altstring+=complement_bases[j];
							}
							else{
								altstring+=bases[j];
							}
						}
					}
					if (NS4string.length()>5){
						NS4string+=",";
					}
					NS4string+=std::to_string(it->second[j]);
					if (isrev){
						j--;
					}
					else{
						j++;
					}
				}
				if (altstring.length()==0){
					altstring=".";
				}
				/*if (i+kmersize+totalbases == 55703){
					std::cout << base << " " << i << std::endl;
				}*/
				char mybase;
				if(isrev){
					mybase=complement(base);
				}
				else{
					mybase=base;
				}

				bool isrepeat=false;
				std::string repeatstring="";
				std::set < std::string >::iterator itb = repeatKmerSet.find(kmer);//check if the kmer is in the hash
				if ( itb != repeatKmerSet.end() ){//if the kmer is in the hash
					isrepeat=true;
					repeatstring=";RR";
				}

				bool isSNP=false;
				if (basecount>1){
					isSNP=true;
				}
				else{
					NS4string="";
				}

				std::string annotation="";

				if ( cait != contigAnnotation.end() ){//if the contig is in the hash
					annotation=features[cait->second[i+kmersize]][8];
				}

				if ((isSNP || not snponly) and (includerepeats || not isrepeat )){
					//std::cout << codons_from_ascii(kmer) << "\t" << base << "\t" << name << "\t" << i+kmersize+1 << "\t" << baseintrep << "\t"<< isSNP << "\t" << isrev << "\t" << isrepeat << std::endl;
					vcffile << name << "\t" << i+kmersize+1 << "\t.\t" << mybase << "\t" << altstring << "\t.\t.\t" << "NS=" << samplesmatching << NS4string << repeatstring << annotation << std::endl;
				}
				
			
			}

		}
		totalbases += sequence.length();
	}
	gzfileStream.close();//close the file
	vcffile.close();

	return 0;
	
}



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
//#include <cctype> //std::tolower
#include "kmers.hpp"
#include "general.hpp"
#include "DNA.hpp"
#include "gzstream.h"


std::unordered_map < std::string, char  > geneticCode = { 
 { "AAA", 'K' }, { "AAG", 'K' },
 { "AAC", 'N' }, { "AAT", 'N' }, 
 { "ACA", 'T' }, { "ACC", 'T' }, { "ACG", 'T' }, { "ACT", 'T' }, 
 { "AGA", 'R' }, { "AGG", 'R' }, { "CGA", 'R' }, { "CGC", 'R' }, { "CGG", 'R' }, { "CGT", 'R' },  
 { "AGC", 'S' }, { "AGT", 'S' }, { "TCA", 'S' }, { "TCC", 'S' }, { "TCG", 'S' }, { "TCT", 'S' }, 
 { "ATA", 'I' }, { "ATC", 'I' }, { "ATT", 'I' },
 { "ATG", 'M' },
 { "CAA", 'Q' }, { "CAG", 'Q' },
 { "CAC", 'H' }, { "CAT", 'H' },
 { "CCA", 'P' }, { "CCC", 'P' }, { "CCG", 'P' }, { "CCT", 'P' },
 { "CTA", 'L' }, { "CTC", 'L' }, { "CTG", 'L' }, { "CTT", 'L' }, { "TTA", 'L' }, { "TTG", 'L' },
 { "GAA", 'E' }, { "GAG", 'E' },
 { "GAC", 'D' }, { "GAT", 'D' },
 { "GCA", 'A' }, { "GCC", 'A' }, { "GCG", 'A' }, { "GCT", 'A' },
 { "GGA", 'G' }, { "GGC", 'G' }, { "GGG", 'G' }, { "GGT", 'G' },
 { "GTA", 'V' }, { "GTC", 'V' }, { "GTG", 'V' }, { "GTT", 'V' },
 { "TAA", '*' }, { "TAG", '*' }, { "TGA", '*' },//Note TGA = W in mycoplasma (translation table 4)
 { "TAC", 'Y' }, { "TAT", 'Y' },
 { "TGC", 'S' }, { "TGT", 'C' },
 { "TGG", 'W' },
 { "TTC", 'F' }, { "TTT", 'F' },
};

std::unordered_map < std::string, char > getGeneticCode(int translationTable){

	std::unordered_map < std::string, char > myGeneticCode=geneticCode;

	if (translationTable==4){
		myGeneticCode.insert(std::make_pair( "TGA", 'W' ));
		
	}
	return myGeneticCode;


}

int findKmersInFasta(const std::vector < std::string > & queryfiles, const std::string & reffile, const bool snponly, const bool includerepeats, const bool includeproduct, const std::string & outputfile)
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
	bool isgff=false;

	int substringlength=(kmersize*2)+1;
	while (gzfileStream.peek() != EOF){

		while (gzfileStream.peek()=='\n'){
			gzfileStream.ignore(256,'\n');//skip the end ofline character
		}

		if (gzfileStream.peek()=='#'){
			std::string gffline;
			while (std::getline(gzfileStream,gffline)){
				if (gffline.length()>1 && gffline[0]=='#' && gffline[1]=='#'){
					if (gffline=="##FASTA"){
						isgff=true;
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
		std::cout << "Reading " << queryfile << std::endl;	
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



	std::unordered_map < std::string, char > geneticcode = getGeneticCode(11);

	//std::cout << repeatKmerSet.size() << " repeat kmers in set" << std::endl;
	
	std::cout << kmerMap.size() << " unique kmers in map" << std::endl;
	
	std::unordered_map < std::string, std::vector < int > > matchMap;

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
	vcffile << "##INFO=<ID=NS5,Number=5,Type=Integer,Description=\"Number Samples With A, C, G, T or N\">" << std::endl;
	if (isgff){
		vcffile << "##INFO=<ID=FT,Number=1,Type=String,Description=\"Feature Type\">" << std::endl;
		vcffile << "##INFO=<ID=FS,Number=1,Type=Character,Description=\"Feature Strand\">" << std::endl;
		vcffile << "##INFO=<ID=BP,Number=1,Type=Integer,Description=\"Position Of Base In Feature\">" << std::endl;
		vcffile << "##INFO=<ID=CP,Number=1,Type=Integer,Description=\"Position Of Base in Codon\">" << std::endl;
		vcffile << "##INFO=<ID=AAP,Number=1,Type=Integer,Description=\"Position Of Amino Acid In Feature\">" << std::endl;
		vcffile << "##INFO=<ID=RAA,Number=1,Type=Character,Description=\"Reference Amino Acid\">" << std::endl;
		vcffile << "##INFO=<ID=AAA,Number=A,Type=Character,Description=\"Alternate Amino Acid(s)\">" << std::endl;
		vcffile << "##INFO=<ID=ID,Number=1,Type=String,Description=\"Feature ID\">" << std::endl;
		vcffile << "##INFO=<ID=LT,Number=1,Type=String,Description=\"Feature Locus Tag\">" << std::endl;
		vcffile << "##INFO=<ID=SI,Number=1,Type=String,Description=\"Feature Systematic ID\">" << std::endl;
		vcffile << "##INFO=<ID=GE,Number=1,Type=String,Description=\"Gene Name\">" << std::endl;
		if (includeproduct){
			vcffile << "##INFO=<ID=PR,Number=1,Type=String,Description=\"Product\">" << std::endl;
		}
	}
	vcffile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;

	int totalbases=0;

	std::cout << "Annotating split kmer matches" << std::endl;
	int matches=0;

	if(openGzFileStream(reffile, gzfileStream)){return 1;}//open the fasta file. May be gzipped.

	std::unordered_map < std::string, std::vector < int > > contigAnnotation;
	std::vector < std::vector < std::string > > features;

	while (gzfileStream.peek() != EOF){

		while (gzfileStream.peek()=='\n'){
			gzfileStream.ignore(256,'\n');//skip the end ofline character
		}

		if (gzfileStream.peek()=='#'){
			std::string gffline;
			while (std::getline(gzfileStream,gffline)){
				if ((gffline.length()==0) ||  ((gffline[0]=='#' && gffline.length()==1) || gffline[0]!='#')){
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

					if (words[2]!="CDS" && words[2]!="tRNA" && words[2]!="rRNA"){
						continue;
					}

					std::unordered_map < std::string, std::vector < int > >::iterator it = contigAnnotation.find(words[0]);//check if the contig is in the hash
					if ( it != contigAnnotation.end() ){//if the contig is in the hash
						if (std::stoi(words[4])>it->second.size()){
							it->second.reserve(std::stoi(words[4]));
						}
						for (int i=std::stoi(words[3])-1; i<std::stoi(words[4]); ++i){
							it->second[i]=features.size();
						}
						//std::cout << it->second.size() << std::endl;
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
			std::cout << features.size() << " features found in " << contigAnnotation.size() << " contigs" << std::endl;
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
		

		//std::cout << "Parsing Sequence" << std::endl;
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
				std::string NS5string=";NS5=";
				std::vector < char > altvector;
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
								altvector.push_back(complement_bases[j]);
							}
							else{
								altstring+=bases[j];
								altvector.push_back(bases[j]);
							}
						}
					}
					if (NS5string.length()>5){
						NS5string+=",";
					}
					NS5string+=std::to_string(it->second[j]);
					if (isrev){
						j--;
					}
					else{
						j++;
					}
				}
				if(it->second[4]>0){//deal with Ns here
					basecount++;
					if (altstring.length()>0){
						altstring+=",";
					}
					altstring+="N";
				}
				if (NS5string.length()>5){
					NS5string+=",";
				}
				NS5string+=std::to_string(it->second[4]);
				samplesmatching+=it->second[4];

				if (altstring.length()==0){
					altstring=".";
				}
				

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
					isrepeat=true;//make it a repeat
					repeatstring=";RR";
				}

				bool isSNP=false;
				if (basecount>1){
					isSNP=true;
				}
				else{
					NS5string="";
				}

				if ((isSNP || not snponly) && (includerepeats || not isrepeat )){
					matches++;
					std::string featureBasePositionString="";
					std::string featureAAPositionString="";
					std::string strandString="";
					std::string codonPositionString="";
					std::string annotationString="";

					std::string annotation="";
					char strand='.';
					int featureStart=-1;
					int featureEnd=-1;
					int featurePhase=0;
					std::string refAA="";
					std::string altAAs="";
					std::string featureType="";
					int featureBasePosition=-1;
					int featureAAPosition=-1;
					int codonPosition=-1;
					

					if ( cait != contigAnnotation.end() ){//if the contig is in the hash
						if (cait->second[i+kmersize]>=0){
							std::string IDString="";
							std::string LTString="";
							std::string SIString="";
							std::string geneString="";
							std::string productString="";
							annotation=features[cait->second[i+kmersize]][8];
							strand=features[cait->second[i+kmersize]][6][0];
							featureStart=stoi(features[cait->second[i+kmersize]][3]);
							featureEnd=stoi(features[cait->second[i+kmersize]][4]);
							featureType=";FT="+features[cait->second[i+kmersize]][2];
							
							if (strand=='+'){
								featureBasePosition=((i+kmersize+1)-featureStart)+1;
							}
							else if (strand=='-'){
								featureBasePosition=(featureEnd-(i+kmersize+1))+1;
							}

							std::stringstream annotationstream;
							annotationstream << annotation;//convert the text to a stringstream
							std::string annotationword;
							while (std::getline(annotationstream, annotationword, ';')){
								std::size_t equalsposition=annotationword.find_first_of("=");
								std::string variable=annotationword.substr(0,equalsposition);
								std::transform(variable.begin(), variable.end(), variable.begin(), ::tolower);
								std::string value=annotationword.substr(equalsposition+1, annotationword.length()-equalsposition);

								value.erase(value.find_last_not_of('"') + 1);
								value.erase(0, value.find_first_not_of('"'));

								if (variable=="gene"){
									geneString=";GE="+value;
								}
								else if (variable=="product"){
									productString=";PR=\""+value+"\"";
								}
								else if (variable=="locus_tag"){
									LTString=";LT="+value;
								}
								else if (variable=="systematic_id"){
									SIString=";SI="+value;
								}
								else if (variable=="id"){
									IDString=";ID="+value;
								}
							}
							annotationString=IDString+LTString+SIString+geneString;
							if (includeproduct){
								annotationString+=productString;
							}
						}
						
					}

					if (featureBasePosition>0){
						featureBasePositionString=";BP="+std::to_string(featureBasePosition);
					}
					if (featureType==";FT=CDS"){
						featurePhase=stoi(features[cait->second[i+kmersize]][7]);
						strandString=";FS=";
						strandString+=strand;
						featureAAPosition=((featureBasePosition+featurePhase-1)/3)+1;
						featureAAPositionString=";AAP="+std::to_string(featureAAPosition);
						codonPosition=(featureBasePosition+featurePhase-1)%3;
						codonPositionString=";CP=";
						codonPositionString+=std::to_string(codonPosition+1);
						std::string featureSequence=sequence.substr(featureStart-1,1+featureEnd-featureStart);
						if (strand=='-'){
							std::transform(featureSequence.begin(),featureSequence.end(),featureSequence.begin(),complement);
							std::reverse(featureSequence.begin(), featureSequence.end());
						}
						std::string refcodon=featureSequence.substr(featureBasePosition-(1+codonPosition),3);
						std::string altcodon;
						char refAAchar=geneticcode[refcodon];
						refAA=";RAA=";
						refAA+=refAAchar;
						std::vector < char > altAAchars;
						for (std::vector < char >::iterator avit=altvector.begin(); avit!=altvector.end(); ++avit){
							altcodon=refcodon;
							altcodon[codonPosition]=*avit;
							altAAchars.push_back(geneticcode[altcodon]);
						}
						if (altAAchars.size()>0){
							altAAs=";AAA=";
							int aacount=0;
							for (std::vector < char >::iterator racit=altAAchars.begin(); racit!=altAAchars.end(); ++racit, ++aacount){
								if (aacount>0){
									altAAs+=",";
								}
								altAAs+=*racit;
							}
						}
					}
				
					vcffile << name << "\t" << i+kmersize+1 << "\t.\t" << mybase << "\t" << altstring << "\t.\t.\t";
					vcffile << "NS=" << samplesmatching << NS5string << repeatstring;
					vcffile << featureType << strandString << featureBasePositionString << codonPositionString << featureAAPositionString << refAA << altAAs << annotationString << std::endl;
				}
				
			
			}

		}
		totalbases += sequence.length();
	}
	gzfileStream.close();//close the file
	vcffile.close();

	std::cout << matches << " split kmer matches annotated" << std::endl << std::endl;
	printDuration(start);

	return 0;
	
}



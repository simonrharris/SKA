#include <iostream> //std::std::cout
#include <fstream> //std::fileStream
#include <sstream> //std::getline
#include <string> //std::string
#include <vector> //std::vector
#include <stdlib.h> //std::strtol
#include <stdexcept> // std::invalid_argument
#include "sk_alleles.hpp"
#include "sk_align.hpp"
#include "sk_compare.hpp"
#include "sk_distance.hpp"
#include "sk_fasta.hpp"
#include "sk_fastq.hpp"
#include "sk_find.hpp"
#include "sk_merge.hpp"
#include "sk_humanise.hpp"
#include "sk_info.hpp"
#include "sk_refalign.hpp"
#include "sk_summary.hpp"
#include "sk_type.hpp"
#include "sk_unique.hpp"
#include "sk_weed.hpp"
#include "general.hpp"

long int getint(char * arg);
float getfloat(char * arg);
int skaHelp(void);
int skaVersion(void);
int allelesHelp(void);
int allelesSubcommand(int argc, char *argv[]);
int alignHelp(void);
int alignSubcommand(int argc, char *argv[]);
int annotateHelp(void);
int annotateSubcommand(int argc, char *argv[]);
int compareHelp(void);
int compareSubcommand(int argc, char *argv[]);
int distanceHelp(void);
int distanceSubcommand(int argc, char *argv[]);
int fastaHelp(void);
int fastaSubcommand(int argc, char *argv[]);
int fastqHelp(void);
int fastqSubcommand(int argc, char *argv[]);
int humaniseHelp(void);
int humaniseSubcommand(int argc, char *argv[]);
int infoHelp(void);
int infoSubcommand(int argc, char *argv[]);
int mapHelp(void);
int mapSubcommand(int argc, char *argv[]);
int mergeHelp(void);
int mergeSubcommand(int argc, char *argv[]);
int summaryHelp(void);
int summarySubcommand(int argc, char *argv[]);
int typeHelp(void);
int typeSubcommand(int argc, char *argv[]);
int uniqueHelp(void);
int uniqueSubcommand(int argc, char *argv[]);
int weedHelp(void);
int weedSubcommand(int argc, char *argv[]);

int MinKmer=3;
int MaxKmer=60;
int MinQual=0;
int MaxQual=60;
int MinCov=1;
int MinFileCov=0;
int maxNamesToPrint=10;


std::string citation = "SKA is currently only available as a preprint, so for now, if you use it, please cite:\n\t\tHarris SR. 2018. SKA: Split Kmer Analysis Toolkit for Bacterial\n\t\tGenomic Epidemiology. bioRxiv 453142 doi: https://doi.org/10.1101/453142";

long int getint(char * arg){
	long int i;
	char * pEnd;i=strtol(arg, &pEnd, 10);
	if(*pEnd != '\0'){
    	throw std::invalid_argument( "received non integer" );
  	}
  	return i;
}

float getfloat(char * arg){
	float f;
	char * pEnd;f=strtof(arg, &pEnd);
	if(*pEnd != '\0'){
    	throw std::invalid_argument( "received non float" );
  	}
  	return f;
}


int skaHelp(void){
	skaVersion();
	std::cout << "Usage:" << std::endl;
	std::cout << "ska <subcommand>" << std::endl << std::endl;
	std::cout << "Subcommands:" << std::endl;
	std::cout << "align\t\tReference-free alignment from a set of split kmer files" << std::endl;
	std::cout << "alleles\t\tCreate a merged split kmer file for all sequenes in one or" << std::endl << "\t\tmore multifasta files" << std::endl;
	std::cout << "annotate\tLocate/annotate split kmers in a reference fasta/gff file" << std::endl;
	std::cout << "compare\t\tPrint comparison statistics for a query split kmer file" << std::endl << "\t\tagainst a set of subject split kmer files" << std::endl;
	std::cout << "distance\tPairwise distance calculation and clustering from split kmer" << std::endl << "\t\tfiles" << std::endl;
	std::cout << "fasta\t\tCreate a split kmer file from fasta file(s)" << std::endl;
	std::cout << "fastq\t\tCreate a split kmer file from fastq file(s)" << std::endl;
	std::cout << "humanise\tPrint kmers from a split kmer file in human readable format" << std::endl;
	std::cout << "info\t\tPrint some information about one or more kmer files" << std::endl;
				          //123456789012345678901234567890123456789012345678901234567890
	std::cout << "map\t\tAlign split kmer file(s) against a reference fasta file" << std::endl;
	std::cout << "merge\t\tMerge split kmer file(s) into one file" << std::endl;
	std::cout << "summary\t\tPrint split kmer file summary statistics" << std::endl;
	std::cout << "type\t\tType split kmer files using a set of allele files" << std::endl;
	std::cout << "unique\t\tOutput kmers unique to a set of split kmer files" << std::endl;
	std::cout << "version\t\tPrint the version and citation for ska" << std::endl;
	std::cout << "weed\t\tWeed kmers from a split kmer file" << std::endl << std::endl;

	               
	return 0;
}

int skaVersion(void){
	std::cout << std::endl << "SKA: Split Kmer Analysis" << std::endl;
	std::cout << "Version: " << versionNumber << std::endl;
	std::cout << "Citation: " << citation << std::endl << std::endl;
	return 0;
}
int allelesHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska fasta [options] <fasta files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-f <file>\tFile of split kmer file names. These will be added to or" << std::endl << "\t\tused as an alternative input to the list provided on the" << std::endl << "\t\tcommand line." << std::endl;
	std::cout << "-k <int>\tSplit Kmer size. The kmer used for searches will be twice" << std::endl << "\t\tthis length, with the variable base in the middle. e.g. a" << std::endl << "\t\tkmer of 15 will search for 31 base matches with the middle" << std::endl << "\t\tbase being allowed to vary. Must be divisible by 3." << std::endl << "\t\t[Default = 15]" << std::endl << std::endl;
	return 0;
}


int allelesSubcommand(int argc, char *argv[]){

	long int kmersize=15;
	std::vector < std::string > args;

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			allelesHelp();
			return 0;
		}
		else if (arg == "-k"){
			i++;
			if (i<argc){
				try {
					kmersize = getint(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting positive integer after -k flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting positive integer after -k flag" << std::endl << std::endl;
				return 1;
			}
			if (kmersize < MinKmer || kmersize > MaxKmer || kmersize % 3 != 0){
				std::cerr << std::endl << "Kmer size must be a positive integer between " << MinKmer << " and " << MaxKmer << " and must be divisible by 3" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl << std::endl;
		allelesHelp();
		return 1;
	}

	skaVersion();

	std::cout << "Creating split " << kmersize << "mers for ";
	if (args.size()>maxNamesToPrint){
		std::cout << args.size() << " files";
	}
	else {
		for (std::vector <std::string>::iterator it = args.begin(); it != args.end(); ++it){
			std::cout << *it << " ";
		}
	}
	std::cout << std::endl << std::endl;

	if (allelesToKmers(args, kmersize)){return 1;}

	return 0;
}



int alignHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska align [options] <split kmer files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << std::endl;
	std::cout << "-k\t\tPrint aligned split kmers to file." << std::endl;
	std::cout << "-o <file>\tPrefix for output files. [Default = reference_free]" << std::endl;
	std::cout << "-p <float>\tMinimum proportion of isolates required to possess a split \n\t\tkmer for that kmer to be included in the alignment. \n\t\t[Default = 0.9]" << std::endl;
	std::cout << "-s <file>\tFile of sample names to include in the alignment." << std::endl;
	std::cout << "-v\t\tOutput variant only alignment. [Default = all sites]" << std::endl << std::endl;
	return 0;
}


int alignSubcommand(int argc, char *argv[]){

	std::string outprefix="reference_free";
	float minproportion=0.9;
	bool variantonly=false;
	bool printkmers=false;
	std::vector < std::string > args;
	std::vector < std::string > sample;

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			alignHelp();
			return 0;
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outprefix = std::string(argv[i]);
			}
			else {
				std::cerr << std::endl << "Expecting string after -o flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-p"){
			i++;
			if (i<argc){
				try {
					minproportion = getfloat(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting float between 0 and 1 after -p flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting float between 0 and 1 after -p flag" << std::endl << std::endl;
				return 1;
			}
			if (minproportion < 0 || minproportion > 1){
				std::cerr << std::endl << "Minimum proportion must be between 0 and 1" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-k"){
			printkmers = true;
		}
		else if (arg == "-s"){
			i++;
			if (i<argc){
				fileToVector(argv[i], sample);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -s flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-v"){
			variantonly = true;
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (variantonly){
		outprefix+="_variants";
	}

	//outprefix+=".aln";
	
	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl;
		alignHelp();
		return 1;
	}

	skaVersion();

	std::cout << "Creating reference free alignment for ";
	if (args.size()>maxNamesToPrint){
		std::cout << args.size() << " files";
	}
	else {	
		for (std::vector <std::string>::iterator it = args.begin(); it != args.end(); ++it){
			std::cout << *it << " ";
		}
	}
	std::cout << std::endl;
	std::cout << "Variable sites present in more than " << minproportion*100 << "% of isolates will be included" << std::endl;
	std::cout << "Output will be written to " << outprefix << std::endl << std::endl;

	if (alignKmers(minproportion, outprefix, args, variantonly, printkmers, sample)){return 1;}

	return 0;
}

int annotateHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska annotate [options] <kmer files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << std::endl;
	std::cout << "-i\t\tInclude kmers in repetitive reference regions." << std::endl;
	std::cout << "-o <file>\tPrefix for output files. [Default = annotation]" << std::endl;
	std::cout << "-p\t\tInclude product in output." << std::endl;
	std::cout << "-r <file>\tReference fasta/gff file name. [Required]" << std::endl;
	std::cout << "-v\t\tOnly output variant sites." << std::endl << std::endl;
	return 0;
}


int annotateSubcommand(int argc, char *argv[]){

	std::string outfile="annotation.vcf";
	std::string reference="";
	std::vector < std::string > args;
	bool snponly=false;
	bool includerepeats=false;
	bool includeproduct=false;

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			annotateHelp();
			return 0;
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-i"){
			includerepeats = true;
		}
		else if (arg == "-p"){
			includeproduct = true;
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = std::string(argv[i])+".vcf";
			}
			else {
				std::cerr << std::endl << "Expecting file name after -o flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-r"){
			i++;
			if (i<argc){
				reference = argv[i];
			}
			else {
				std::cerr << std::endl << "Expecting file name after -r flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-v"){
			snponly = true;
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl;
		annotateHelp();
		return 1;
	}
	if (reference==""){
		std::cerr << std::endl << "Reference sequence is required" << std::endl;
		annotateHelp();
		return 1;
	}

	skaVersion();

	std::cout << "Annotating split kmer matches in " << reference << " from ";
	if (args.size()>maxNamesToPrint){
		std::cout << args.size() << " files";
	}
	else {	
		for (std::vector <std::string>::iterator it = args.begin(); it != args.end(); ++it){
			std::cout << *it << " ";
		}
	}
	std::cout << std::endl;

	if (snponly){
		std::cout << "Only variant sites will be included" << std::endl;
	}
	if (includerepeats){
		std::cout << "Matches in repeats will be included" << std::endl;
	}
	else {
		std::cout << "Matches in repeats will not be included" << std::endl;
	}
	if (includeproduct){
		std::cout << "Product information for CDSs will be printed" << std::endl;
	}
	else{
		std::cout << "Product information for CDSs will not be printed" << std::endl;
	}

	std::cout << "Output will be written to " << outfile << std::endl << std::endl;


	if (findKmersInFasta(args, reference, snponly, includerepeats, includeproduct, outfile)){return 1;}

	return 0;
}

int compareHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska compare [options] <subject split kmer files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << std::endl;
	std::cout << "-q <file>\tQuery split kmer file." << std::endl << std::endl;
	return 0;
}


int compareSubcommand(int argc, char *argv[]){

	std::vector < std::string > args;
	std::string queryfile="";

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			compareHelp();
			return 0;
		}
		else if (arg == "-q"){
			i++;
			if (i<argc){
				queryfile = argv[i];
			}
			else {
				std::cerr << std::endl << "Expecting file name after -q flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl;
		compareHelp();
		return 1;
	}

	if (queryfile==""){
		std::cerr << std::endl << "Query kmer file is required" << std::endl;
		compareHelp();
		return 1;
	}

	if (compareKmerFiles(queryfile, args)){return 1;}

	return 0;
}

int distanceHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska distance [options] <split kmer files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-c \t\tDo not print clusters files." << std::endl;
	std::cout << "-d \t\tDo not print distances file." << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << std::endl;
	std::cout << "-i <float>\tIdentity cutoff for defining clusters. Isolates will be \n\t\tclustered if they share at least this proportion of the \n\t\tsplit kmers in the file with fewer kmers and pass the SNP \n\t\tcutoff. [Default = 0.9]" << std::endl;
	std::cout << "-o <file>\tPrefix for output files. [Default = distances]" << std::endl;
	std::cout << "-s <int>\tSNP cutoff for defining clusters. Isolates will be clustered \n\t\tif they are separated by fewer than this number of SNPs and \n\t\tpass the identity cutoff. [Default = 20]" << std::endl;
	std::cout << "-S \t\tInclude singletons in dot file" << std::endl << std::endl;
	return 0;
}


int distanceSubcommand(int argc, char *argv[]){

	bool distancefile=true;
	bool clusterfile=true;
	bool includesingletons=false;
	std::string prefix="distances";
	float minid=0.9;
	int maxsnps=20;
	std::vector < std::string > args;

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			distanceHelp();
			return 0;
		}
		else if (arg=="-c"){
			clusterfile=false;
		}
		else if (arg == "-d"){
			distancefile=false;
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-i"){
			i++;
			if (i<argc){
				try {
					minid = getfloat(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting float between 0 and 1 after -i flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting float between 0 and 1 after -i flag" << std::endl << std::endl;
				return 1;
			}
			if (minid < 0 || minid > 1){
				std::cerr << std::endl << "Minimum identity must be between 0 and 1" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				prefix = argv[i];
			}
			else {
				std::cerr << std::endl << "Expecting string after -o flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-s"){
			i++;
			if (i<argc){
				try {
					maxsnps = getint(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting positive integer after -s flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting positive integer after -s flag" << std::endl << std::endl;
				return 1;
			}
			if (maxsnps < 0){
				std::cerr << std::endl << "Maximum number of SNPs must be 0 or greater" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-S"){
			includesingletons=true;
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}
	
	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl;
		distanceHelp();
		return 1;
	}

	if (distancefile==false && clusterfile==false){
		std::cerr << std::endl << "At least one output file is required. i.e. you cannot use the -c and -d flags together." << std::endl << std::endl;
		return 1;
	}
	if (prefix=="false"){
		std::cerr << std::endl << "An output prefix (-p) is required." << std::endl << std::endl;
		return 1;
	}

	skaVersion();

	std::cout << "Calculating pairwise distances for ";
	if (args.size()>maxNamesToPrint){
		std::cout << args.size() << " files";
	}
	else {
		for (std::vector < std::string >::iterator it = args.begin(); it != args.end(); ++it){
			std::cout << *it << " ";
		}
	}
	std::cout << std::endl;


	std::cout << "Clusters will be created for isolates that are within " << maxsnps << " SNPs of one another and share at least " << minid*100 << "% of the split kmers of the isolate with fewer kmers" << std::endl;

	if (distancefile){
		std::cout << "Distances will be written to " << prefix << ".distances.tsv" << std::endl;
	}
	if (clusterfile){
		std::cout << "Clusters will be written to files with the prefix " << prefix << std::endl;
	}
	std::cout << std::endl;

	if (kmerDistance(prefix, distancefile, clusterfile, args, maxsnps, minid, includesingletons)){return 1;}

	return 0;
}

int fastaHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska fasta [options] <fasta files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-c\t\tTreat all contigs as circular." << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << std::endl;
	std::cout << "-k <int>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\t[Default = 15]" << std::endl;
	std::cout << "-o <file>\tOutput file prefix. [Default = fasta]" << std::endl << std::endl;
	return 0;
}


int fastaSubcommand(int argc, char *argv[]){

	std::string outfile="fasta.skf";
	long int kmersize=15;
	bool circular=false;
	std::vector < std::string > args;

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			fastaHelp();
			return 0;
		}
		else if (arg=="-c"){
			circular=true;
		}
		else if (arg == "-k"){
			i++;
			if (i<argc){
				try {
					kmersize = getint(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting positive integer after -k flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting positive integer after -k flag" << std::endl << std::endl;
				return 1;
			}
			if (kmersize < MinKmer || kmersize > MaxKmer || kmersize % 3 != 0){
				std::cerr << std::endl << "Kmer size must be a positive integer between " << MinKmer << " and " << MaxKmer << " and must be divisible by 3" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = std::string(argv[i])+".skf";
			}
			else {
				std::cerr << std::endl << "Expecting file name after -o flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl;
		fastaHelp();
		return 1;
	}

	skaVersion();

	std::cout << "Creating split " << kmersize << "mers for ";
	if (args.size()>maxNamesToPrint){
		std::cout << args.size() << " files";
	}
	else {
		for (std::vector < std::string >::iterator it = args.begin(); it != args.end(); ++it){
			std::cout << *it << " ";
		}
	}
	std::cout << std::endl;
	if (circular){
		std::cout << "Contigs will be treated as circular" << std::endl;
	}
	else {
		std::cout << "Contigs will be treated as linear" << std::endl;
	}
	std::cout << "Output will be written to " << outfile << std::endl;

	std::cout << std::endl;

	if (fastaToKmers(args, outfile, kmersize, circular)){return 1;}

	return 0;
}



int fastqHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska fastq [options] <fastq files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-a\t\tPrint allele frequencies of split kmers to file [Default = false]" << std::endl;
	std::cout << "-c <int>\tCoverage cutoff. Kmers with coverage below this value will \n\t\tbe discarded. [Default = 4]" << std::endl;
	std::cout << "-C <int>\tFile coverage cutoff. Kmers with coverage below this value \n\t\tin any of the fastq files will be discarded. [Default = 2]" << std::endl;
	std::cout << "-k <int>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\t[Default = 15]" << std::endl;
	std::cout << "-m <float>\tMinimum allowable minor allele frequency. Kmer alleles below \n\t\tthis frequency will be discarded. [Default = 0.2]" << std::endl;
	std::cout << "-o <file>\tOutput prefix. [Default = fastq]" << std::endl;
	std::cout << "-q <int>\tQuality filter for fastq files. No kmers will be created \n\t\tfrom sequence including quality scores below this cutoff. \n\t\t[Default = 20]" << std::endl << std::endl;
	return 0;
}


int fastqSubcommand(int argc, char *argv[]){

	std::string outfile="fastq";
	long int covcutoff=4;
	long int filecovcutoff=2;
	long int kmersize=15;
	long int minquality=20;
	bool printalleles=false;
	float minmaf=0.2;
	std::vector < std::string > args;

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			fastqHelp();
			return 0;
		}
		else if (arg=="-a"){
			printalleles=true;
		}
		else if (arg == "-c"){
			i++;
			if (i<argc){
				try {
					covcutoff = getint(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting positive integer after -c flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting positive integer after -c flag" << std::endl << std::endl;
				return 1;
			}
			if (covcutoff < MinCov){
				std::cerr << std::endl << "Coverage cutoff must be a positive integer greater than or equal to " << MinCov << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-C"){
			i++;
			if (i<argc){
				try {
					filecovcutoff = getint(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting positive integer after -C flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting positive integer after -C flag" << std::endl << std::endl;
				return 1;
			}
			if (filecovcutoff < MinFileCov){
				std::cerr << std::endl << "File coverage cutoff must be a positive integer greater than or equal to " << MinFileCov << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-k"){
			i++;
			if (i<argc){
				try {
					kmersize = getint(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting positive integer after -k flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting positive integer after -k flag" << std::endl << std::endl;
				return 1;
			}
			if (kmersize < MinKmer || kmersize > MaxKmer || kmersize % 3 != 0){
				std::cerr << std::endl << "Kmer size must be a positive integer between " << MinKmer << " and " << MaxKmer << " and must be divisible by 3" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-m"){
			i++;
			if (i<argc){
				try {
					minmaf = getfloat(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting float between 0 and 1 after -m flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting float between 0 and 1 after -m flag" << std::endl << std::endl;
				return 1;
			}
			if (minmaf < 0 || minmaf > 1){
				std::cerr << std::endl << "Minimum minor allele frequency must be between 0 and 1" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = std::string(argv[i]);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -o flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-q"){
			i++;
			if (i<argc){
				try {
					minquality = getint(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting positive integer after -q flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting positive integer after -q flag" << std::endl << std::endl;
				return 1;
			}
			if (minquality < MinQual || minquality > MaxQual){
				std::cerr << std::endl << "Minimum quality cutoff must be a positive integer between " << MinQual << " and " << MaxQual << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl;
		fastqHelp();
		return 1;
	}

	skaVersion();

	std::cout << "Creating split " << kmersize << "mers for ";
	if (args.size()>maxNamesToPrint){
		std::cout << args.size() << " files";
	}
	else {
		for (std::vector < std::string >::iterator it = args.begin(); it != args.end(); ++it){
			std::cout << *it << " ";
		}
	}
	std::cout << std::endl;
	std::cout << "Filters:" << std::endl;
	std::cout << "\tCoverage cutoff = " << covcutoff << std::endl;
	std::cout << "\tPer file coverage cutoff = " << filecovcutoff << std::endl;
	std::cout << "\tMinimum minor allele frequency = " << minmaf << std::endl;
	std::cout << "\tMinimum base quality = " << minquality << std::endl;
	std::cout << "Output will be written to " << outfile << ".skf" << std::endl;
	if (printalleles){
		std::cout << "Allele frequencies will be written to " << outfile << "_alleles.tsv" << std::endl;
	}

	std::cout << std::endl;

	if (fastqToKmers(args, outfile, kmersize, minquality, filecovcutoff, covcutoff, minmaf, printalleles)){return 1;}

	return 0;
}


int humaniseHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska humanise [options]" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-i <file>\tInput split kmer file." << std::endl;
	std::cout << "-o <file>\tOutput prefix. [Default = humanised_kmers]" << std::endl;
	std::cout << "-t\t\tPrint tabulated output." << std::endl << std::endl;
	return 0;
}


int humaniseSubcommand(int argc, char *argv[]){

	std::string inputfile="";
	std::string outfile="humanised_kmers.tsv";

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			infoHelp();
			return 0;
		}
		else if (arg == "-i"){
			i++;
			if (i<argc){
				inputfile = argv[i];
			}
			else {
				std::cerr << std::endl << "Expecting file name after -i flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = std::string(argv[i])+".tsv";
			}
			else {
				std::cerr << std::endl << "Expecting file name after -o flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
			return 1;
		}
		else {
			std::cerr << std::endl << "Unrecognised argument " << arg << std::endl << std::endl;
			return 1;
		}
	}


	if (inputfile==""){
		std::cerr << std::endl << "Input kmer file is required" << std::endl;
		humaniseHelp();
		return 1;
	}

	if (humaniseKmers(inputfile, outfile)){return 1;}

	return 0;
}


int infoHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska info [options] <split kmer files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << std::endl;
	std::cout << "-t\t\tPrint tabulated output." << std::endl << std::endl;
	return 0;
}


int infoSubcommand(int argc, char *argv[]){

	std::vector < std::string > args;
	bool tabulated=false;

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			infoHelp();
			return 0;
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-t"){
			tabulated=true;
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl;
		infoHelp();
		return 1;
	}

	if (getKmerFileInfo(args, tabulated)){return 1;}

	return 0;
}

int mapHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska map [options] <split kmer files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-a <file>\tMap all bases of kmers (Default = just map middle base)." << std::endl;
	std::cout << "-c\t\tTreat all reference contigs as circular." << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << std::endl;
	std::cout << "-k <int>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\tMust be the same value used to create the kmer files. \n\t\t[Default = 15]" << std::endl;
	std::cout << "-i\t\tInclude reference sequence in alignment." << std::endl;
	std::cout << "-m\t\tMap bases to repeats rather than making them N." << std::endl;
	std::cout << "-o <file>\tOutput file prefix. [Default = mappedkmers]" << std::endl;
	std::cout << "-r <file>\tReference fasta file name. [Required]" << std::endl;
	std::cout << "-s <file>\tFile of sample names to include in the alignment." << std::endl;
	std::cout << "-v\t\tOutput variant only alignment. [Default = all sites]" << std::endl << std::endl;
	return 0;
}


int mapSubcommand(int argc, char *argv[]){

	std::string outprefix = "mappedkmers";
	std::string reference = "";
	long int kmersize=15;
	bool includeref=false;
	bool maprepeats=false;
	bool fillall=false;
	bool variantonly=false;
	bool circular=false;
	std::vector < std::string > sample;
	std::vector < std::string > args;

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			mapHelp();
			return 0;
		}
		else if (arg=="-c"){
			circular=true;
		}
		else if (arg == "-k"){
			i++;
			if (i<argc){
				try {
					kmersize = getint(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting positive integer after -k flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting positive integer after -k flag" << std::endl << std::endl;
				return 1;
			}
			if (kmersize < MinKmer || kmersize > MaxKmer || kmersize % 3 != 0){
				std::cerr << std::endl << "Kmer size must be a positive integer between " << MinKmer << " and " << MaxKmer << " and must be divisible by 3" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-r"){
			i++;
			if (i<argc){
				reference = argv[i];
			}
			else {
				std::cerr << std::endl << "Expecting file name after -r flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-i"){
			includeref = true;
		}
		else if (arg == "-m"){
			maprepeats = true;
		}
		else if (arg == "-a"){
			fillall = true;
		}
		else if (arg == "-v"){
			variantonly = true;
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outprefix = std::string(argv[i]);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -o flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-s"){
			i++;
			if (i<argc){
				fileToVector(argv[i], sample);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -s flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl;
		mapHelp();
		return 1;
	}

	if (variantonly){
		outprefix+="_variants";
	}

	outprefix+=".aln";

	if (reference==""){
		std::cerr << std::endl << "Reference file is required" << std::endl;
		mapHelp();
		return 1;
	}

	skaVersion();

	std::cout << "Aligning ";
	if (args.size()>maxNamesToPrint){
		std::cout << args.size() << " files";
	}
	else {
		for (std::vector < std::string >::iterator it = args.begin(); it != args.end(); ++it){
			std::cout << *it << " ";
		}
	}
	std::cout << std::endl;
	std::cout << "Using split kmer size of " << kmersize << std::endl;
	std::cout << "Using " << reference << " as reference" << std::endl;
	if (includeref){
		std::cout << "Reference will be included in the alignment" << std::endl;
	}
	else{
		std::cout << "Reference will not be included in the alignment" << std::endl;
	}
	if (maprepeats){
		std::cout << "Repeats in the reference will be mapped rather than left as gaps" << std::endl;
	}
	else{
		std::cout << "Repeats in the reference will be output as gaps" << std::endl;
	}
	if (fillall){
		std::cout << "All bases in kmers will be mapped as lower case" << std::endl;
	}
	if (circular){
		std::cout << "Reference contigs will be treated as circular" << std::endl;
	}
	else {
		std::cout << "Reference contigs will be treated as linear" << std::endl;
	}
	std::cout << "Output will be written to " << outprefix << std::endl << std::endl;

	if (alignKmersToReference(reference, outprefix, args, kmersize, includeref, maprepeats, fillall, variantonly, sample, circular)){return 1;}

	return 0;
}


int mergeHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska merge [options] <subject split kmer files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << std::endl;
	std::cout << "-o <file>\tOutput file prefix. [Default = merged]" << std::endl;
	std::cout << "-s <file>\tFile of sample names to include in the merged file. Can \n\t\tbe used to extract samples from a merge" << std::endl << std::endl;
	return 0;
}


int mergeSubcommand(int argc, char *argv[]){

	std::vector < std::string > args;
	std::vector < std::string > sample;
	std::string outfile="merged.skf";

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			mergeHelp();
			return 0;
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = std::string(argv[i])+".skf";
			}
			else {
				std::cerr << std::endl << "Expecting file name after -o flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-s"){
			i++;
			if (i<argc){
				fileToVector(argv[i], sample);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -s flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl;
		mergeHelp();
		return 1;
	}

	if (outfile==""){
		std::cerr << std::endl << "Output file name is required" << std::endl;
		mergeHelp();
		return 1;
	}

	skaVersion();

	std::cout << "Merging ";
	if (args.size()>maxNamesToPrint){
		std::cout << args.size() << " files";
	}
	else {
		for (std::vector < std::string >::iterator it = args.begin(); it != args.end(); ++it){
			std::cout << *it << " ";
		}
	}
	std::cout << std::endl;
	std::cout << "Output will be written to " << outfile << std::endl << std::endl;

	if (mergeKmerFiles(outfile, args, sample)){return 1;}

	return 0;
}


int summaryHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska summary [options] <split kmer files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl << std::endl;
	return 0;
}


int summarySubcommand(int argc, char *argv[]){

	std::vector < std::string > args;

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			summaryHelp();
			return 0;
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl;
		summaryHelp();
		return 1;
	}

	if (summariseKmerFiles(args)){return 1;}

	return 0;
}


int typeHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska type [options] <locus fasta files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-f <file>\tFile of locus fasta file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << std::endl;
	std::cout << "-p <file>\ttab file containing profile information." << std::endl;
	std::cout << "-q <file>\tQuery split kmer file. This can be a single kmer file." << std::endl << std::endl;
	return 0;
}


int typeSubcommand(int argc, char *argv[]){

	std::vector < std::string > args;
	std::string queryfile="";
	std::string profilefile="";

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			typeHelp();
			return 0;
		}
		else if (arg == "-q"){
			i++;
			if (i<argc){
				queryfile = argv[i];
			}
			else {
				std::cerr << std::endl << "Expecting file name after -q flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-p"){
			i++;
			if (i<argc){
				profilefile = argv[i];
			}
			else {
				std::cerr << std::endl << "Expecting file name after -p flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl;
		typeHelp();
		return 1;
	}

	if (queryfile==""){
		std::cerr << std::endl << "Query kmer file is required" << std::endl;
		typeHelp();
		return 1;
	}

	if (typeKmerFile(queryfile, profilefile, args)){return 1;}

	return 0;
}


int uniqueHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska unique [options]" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-i <file>\tFile of ingroup sample names. Unique kmers found \n\t\tin these files will be retained." << std::endl;
	std::cout << "-n\t\tAllow Ns as in unique split kmers." << std::endl;
	std::cout << "-o <file>\tOutput file prefix. [Default = unique]" << std::endl;
	std::cout << "-p <float>\tMinimum proportion of ingroup isolates required to possess a \n\t\tsplit kmer for that kmer to be retained. [Default = 0.9]" << std::endl << std::endl;
					   
	return 0;
}


int uniqueSubcommand(int argc, char *argv[]){

	std::string outfile="unique.skf";
	std::vector < std::string > ingroup;
	float minproportion=0.9;
	std::vector < std::string > args;
	bool allowNs = false;

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			uniqueHelp();
			return 0;
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = std::string(argv[i])+".skf";
			}
			else {
				std::cerr << std::endl << "Expecting file name after -o flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-i"){
			i++;
			if (i<argc){
				fileToVector(argv[i], ingroup);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -i flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-n"){
			allowNs = true;
		}
		else if (arg == "-p"){
			i++;
			if (i<argc){
				try {
					minproportion = getfloat(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting float between 0 and 1 after -p flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting float between 0 and 1 after -p flag" << std::endl << std::endl;
				return 1;
			}
			if (minproportion < 0 || minproportion > 1){
				std::cerr << std::endl << "Minimum proportion must be between 0 and 1" << std::endl << std::endl;
				return 1;
			}
		}
		else {
			args.push_back(arg);
		}
	}

	if (ingroup.size()<1){
		std::cerr << std::endl << "At least one ingroup kmer file must be specified" << std::endl;
		uniqueHelp();
		return 1;
	}

	if (args.size()<1){
		std::cerr << std::endl << "At least one kmer file must be specified" << std::endl;
		uniqueHelp();
		return 1;
	}

	skaVersion();

	std::cout << "Finding unique kmers in at least " << minproportion*100 << "% of ";
	for (std::vector < std::string >::iterator it = ingroup.begin(); it != ingroup.end(); ++it){
		std::cout << *it << " ";
	}

	std::cout << std::endl;

	std::cout << "Output will be written to " << outfile << std::endl;

	std::cout << std::endl;

	if (uniqueKmers(ingroup, args, minproportion, outfile, allowNs)){return 1;}

	return 0;
}

int weedHelp(void){
	std::cout << std::endl << "Usage:" << std::endl;
	std::cout << "ska weed [options] <split kmer files>" << std::endl << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << std::endl;
	std::cout << "-h\t\tPrint this help." << std::endl;
	std::cout << "-i <file>\tName of kmer file containing kmers to be weeded. [Required]" << std::endl;
	std::cout << "-m <int>\tMinimum number of samples required to possess a split\n\t\tkmer for that kmer to be retained. [Default = 0]" << std::endl;
	std::cout << "-M <int>\tMaximum number of samples required to possess a split\n\t\tkmer for that kmer to be retained. 0 = No maximum. [Default = 0]" << std::endl;
	std::cout << "-p <float>\tMinimum proportion of samples required to possess a split\n\t\tkmer for that kmer to be retained. [Default = 0.0]" << std::endl;
	std::cout << "-P <float>\tMaximum proportion of samples required to possess a split\n\t\tkmer for that kmer to be retained. [Default = 1.0]" << std::endl;
	                        //123456789012345678901234567890123456789012345678901234567890
	return 0;
}


int weedSubcommand(int argc, char *argv[]){

	std::string infile="";
	std::string kmerfile="";
	float minproportion=0.0;
	float maxproportion=1.0;
	int minsamples=0;
	int maxsamples=0;
	std::vector < std::string > args;

	for (int i=2; i<argc; ++i){

		std::string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			weedHelp();
			return 0;
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				std::cerr << std::endl << "Expecting file name after -f flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-i"){
			i++;
			if (i<argc){
				infile = argv[i];
			}
			else {
				std::cerr << std::endl << "Expecting file name after -i flag" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-m"){
			i++;
			if (i<argc){
				try {
					minsamples = getint(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting positive int after -m flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting positive int after -m flag" << std::endl << std::endl;
				return 1;
			}
			if (minsamples < 0){
				std::cerr << std::endl << "Minimum number of samples must be 0 or greater" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-M"){
			i++;
			if (i<argc){
				try {
					maxsamples = getint(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting positive int after -m flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting positive int after -m flag" << std::endl << std::endl;
				return 1;
			}
			if (maxsamples < 0){
				std::cerr << std::endl << "Maximum number of samples must be 0 or greater" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-p"){
			i++;
			if (i<argc){
				try {
					minproportion = getfloat(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting float between 0 and 1 after -p flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting float between 0 and 1 after -p flag" << std::endl << std::endl;
				return 1;
			}
			if (minproportion < 0 || minproportion > 1){
				std::cerr << std::endl << "Minimum proportion must be between 0 and 1" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg == "-P"){
			i++;
			if (i<argc){
				try {
					maxproportion = getfloat(argv[i]);
				}
				catch (const std::invalid_argument& e) {
					std::cerr << std::endl << "Expecting float between 0 and 1 after -P flag" << std::endl << std::endl;
					return 1;
				}
			}
			else {
				std::cerr << std::endl << "Expecting float between 0 and 1 after -P flag" << std::endl << std::endl;
				return 1;
			}
			if (maxproportion < 0 || maxproportion > 1){
				std::cerr << std::endl << "Maximum proportion must be between 0 and 1" << std::endl << std::endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag " << arg << std::endl << std::endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		std::cerr << std::endl << "Too few arguments" << std::endl;
		weedHelp();
		return 1;
	}

	if (infile=="" && minproportion==0 && minsamples==0 && maxproportion==1 && maxproportion==0){
		std::cerr << "Input split kmer file is required if no sample numbers or proportions are set" << std::endl;
		weedHelp();
		return 1;
	}

	if (minproportion>=maxproportion){
		std::cerr << "Maxproportion must be greater than minproportion" << std::endl;
		weedHelp();
		return 1;
	}
	if (minsamples>=maxsamples && maxsamples!=0){
		std::cerr << "Maxsamples must be greater than minsamples" << std::endl;
		weedHelp();
		return 1;
	}

	skaVersion();

	std::cout << "Weeding kmers in " << infile << std::endl;
	std::cout << "From ";
	if (args.size()>maxNamesToPrint){
		std::cout << args.size() << " files";
	}
	else {
		for (std::vector < std::string >::iterator it = args.begin(); it != args.end(); ++it){
			std::cout << *it << " ";
		}
	}

	std::cout << std::endl;

	if (weedKmers(args, infile, minproportion, maxproportion, minsamples, maxsamples)){return 1;}

	return 0;
}


int main(int argc, char *argv[])
{
	
	if (argc<2){
		skaHelp();
		return 0;
	}

	std::string subcommand=(argv[1]);

	if (subcommand == "alleles"){
			if(allelesSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "align"){
			if(alignSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "compare"){
			if(compareSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "distance"){
			if(distanceSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "fasta"){
			if(fastaSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "fastq"){
			if(fastqSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "annotate"){
			if(annotateSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "humanise"){
			if(humaniseSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "info"){
			if(infoSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "map"){
			if(mapSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "merge"){
			if(mergeSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "summary"){
			if(summarySubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "type"){
			if(typeSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "unique"){
			if(uniqueSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "weed"){
			if(weedSubcommand(argc, argv)){return 1;}
		}
	else if (subcommand == "-v" || subcommand == "--version" || subcommand == "version"){
			if(skaVersion()){return 1;}
		}
	else if (subcommand == "-h" || subcommand == "--help"){
			if(skaHelp()){return 1;}
		}
	else if (subcommand[0]=='-'){
			std::cerr << std::endl << "Unrecognised flag: " << subcommand << std::endl << std::endl;
			return 1;
		}
	else {
		std::cerr << std::endl << "Unrecognised subcommand: " << subcommand << std::endl;
		skaHelp();
		return 1;
	}
	
	return 0;

}
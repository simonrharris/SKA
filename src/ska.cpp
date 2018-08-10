#include <iostream> //std::cout
#include <fstream> //std::fileStream
#include <sstream> //std::getline
#include <string> //std::string
#include <vector> //std::vector
#include <stdlib.h> //std::strtol
#include <stdexcept>      // std::invalid_argument
#include "sk_alleles.hpp"
#include "sk_align.hpp"
#include "sk_compare.hpp"
#include "sk_distance.hpp"
#include "sk_fasta.hpp"
#include "sk_fastq.hpp"
#include "sk_merge.hpp"
#include "sk_info.hpp"
#include "sk_refalign.hpp"
#include "sk_summary.hpp"
#include "sk_type.hpp"
#include "sk_unique.hpp"
#include "sk_weed.hpp"
#include "general.hpp"

using namespace std;

long int getint(char * arg);
float getfloat(char * arg);
int skaHelp(void);
int skaVersion(void);
int allelesHelp(void);
int allelesSubcommand(int argc, char *argv[]);
int alignHelp(void);
int alignSubcommand(int argc, char *argv[]);
int compareHelp(void);
int compareSubcommand(int argc, char *argv[]);
int distanceHelp(void);
int distanceSubcommand(int argc, char *argv[]);
int fastaHelp(void);
int fastaSubcommand(int argc, char *argv[]);
int fastqHelp(void);
int fastqSubcommand(int argc, char *argv[]);
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
int MaxKmer=30;
int MinQual=0;
int MaxQual=60;
int MinCov=0;
int maxNamesToPrint=10;

string versionNumber = "0.1";
string citation = "TBA";

long int getint(char * arg){
	long int i;
	char * pEnd;i=strtol(arg, &pEnd, 10);
	if(*pEnd != '\0'){
    	throw invalid_argument( "received non integer" );
  	}
  	return i;
}

float getfloat(char * arg){
	float f;
	char * pEnd;f=strtof(arg, &pEnd);
	if(*pEnd != '\0'){
    	throw invalid_argument( "received non float" );
  	}
  	return f;
}


int skaHelp(void){
	skaVersion();
	cout << "Usage:" << endl;
	cout << "ska <subcommand>" << endl << endl;
	cout << "Subcommands:" << endl;
	cout << "align\t\tReference-free alignment from a set of split kmer files" << endl;
	cout << "alleles\t\tCreate a merged split kmer file for all sequenes in one or" << endl << "\t\tmore multifasta files" << endl;
	cout << "compare\t\tPrint comparison statistics for a query split kmer file" << endl << "\t\tagainst a set of subject split kmer files" << endl;
	cout << "distance\tPairwise distance calculation and clustering from split kmer" << endl << "\t\tfiles" << endl;
	cout << "fasta\t\tCreate a split kmer file from fasta file(s)" << endl;
	cout << "fastq\t\tCreate a split kmer file from fastq file(s)" << endl;
	cout << "info\t\tPrint some information about one or more kmer or kmerge" << endl << "\t\tfiles" << endl;
				   //123456789012345678901234567890123456789012345678901234567890
	cout << "map\t\tAlign split kmer file(s) against a reference fasta file" << endl;
	cout << "summary\t\tPrint split kmer file summary statistics" << endl;
	cout << "type\t\tType split kmer files using a set of allele files" << endl;
	cout << "unique\t\tOutput kmers unique to a set of split kmer files" << endl;
	cout << "version\t\tPrint the version and citation for ska" << endl;
	cout << "weed\t\tWeed kmers from a split kmer file" << endl << endl;

	               
	return 0;
}

int skaVersion(void){
	cout << endl << "SKA: Split Kmer Analysis" << endl;
	cout << "Version: " << versionNumber << endl;
	cout << "Citation: " << citation << endl << endl;
	return 0;
}
int allelesHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska fasta [options] <fasta files>" << endl << endl;
	cout << "Options:" << endl;
	cout << "-h\t\tPrint this help." << endl;
	cout << "-f <file>\tFile of split kmer file names. These will be added to or" << endl << "\t\tused as an alternative input to the list provided on the" << endl << "\t\tcommand line." << endl;
	cout << "-k <int>\tSplit Kmer size. The kmer used for searches will be twice" << endl << "\t\tthis length, with the variable base in the middle. e.g. a" << endl << "\t\tkmer of 15 will search for 31 base matches with the middle" << endl << "\t\tbase being allowed to vary. Must be divisible by 3." << endl << "\t\t[Default = 15]" << endl << endl;
	return 0;
}


int allelesSubcommand(int argc, char *argv[]){

	long int kmersize=15;
	vector<string> args;

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

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
				catch (const invalid_argument& e) {
					cerr << endl << "Expecting positive integer after -k flag" << endl << endl;
					return 1;
				}
			}
			else {
				cerr << endl << "Expecting positive integer after -k flag" << endl << endl;
				return 1;
			}
			if (kmersize < MinKmer || kmersize > MaxKmer || kmersize % 3 != 0){
				cerr << endl << "Kmer size must be a positive integer between " << MinKmer << " and " << MaxKmer << " and must be divisible by 3" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cerr << endl << "Expecting file name after -f flag" << endl << endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			cerr << endl << "Unrecognised flag " << arg << endl << endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cerr << endl << "Too few arguments" << endl << endl;
		allelesHelp();
		return 1;
	}

	skaVersion();

	cout << "Creating split " << kmersize << "mers for ";
	if (args.size()>maxNamesToPrint){
		cout << args.size() << " files";
	}
	else {
		for (auto it = args.begin(); it != args.end(); ++it){
			cout << *it << " ";
		}
	}
	cout << endl << endl;

	allelesToKmers(args, kmersize);

	return 0;
}



int alignHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska align [options] <split kmer files>" << endl << endl;
	cout << "Options:" << endl;
	cout << "-h\t\tPrint this help." << endl;
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << endl;
	cout << "-o <file>\tPrefix for output files. [Default = reference_free]" << endl;
	cout << "-p <float>\tMinimum proportion of isolates required to possess a split \n\t\tkmer for that kmer to be included in the alignment. \n\t\t[Default = 0.9]" << endl;
	cout << "-s <file>\tFile of sample names to include in the alignment." << endl;
	cout << "-v\tOutput variant only alignment. [Default = all sites]" << endl << endl;
	return 0;
}


int alignSubcommand(int argc, char *argv[]){

	string outprefix="reference_free";
	float minproportion=0.9;
	bool variantonly=false;
	vector <string> args;
	vector <string> sample;

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			alignHelp();
			return 0;
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outprefix = string(argv[i]);
			}
			else {
				cerr << endl << "Expecting string after -o flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-p"){
			i++;
			if (i<argc){
				try {
					minproportion = getfloat(argv[i]);
				}
				catch (const invalid_argument& e) {
					cerr << endl << "Expecting float between 0 and 1 after -p flag" << endl << endl;
					return 1;
				}
			}
			else {
				cerr << endl << "Expecting float between 0 and 1 after -p flag" << endl << endl;
				return 1;
			}
			if (minproportion < 0 || minproportion > 1){
				cerr << endl << "Minimum proportion must be between 0 and 1" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cerr << endl << "Expecting file name after -f flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-s"){
			i++;
			if (i<argc){
				fileToVector(argv[i], sample);
			}
			else {
				cerr << endl << "Expecting file name after -s flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-v"){
			variantonly = true;
		}
		else if (arg[0]=='-'){
			cerr << endl << "Unrecognised flag " << arg << endl << endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (variantonly){
		outprefix+="_variants";
	}

	outprefix+=".aln";
	
	if (args.size()<1){
		cerr << endl << "Too few arguments" << endl;
		alignHelp();
		return 1;
	}

	skaVersion();

	cout << "Creating reference free alignment for ";
	if (args.size()>maxNamesToPrint){
		cout << args.size() << " files";
	}
	else {	
		for (auto it = args.begin(); it != args.end(); ++it){
			cout << *it << " ";
		}
	}
	cout << endl;
	cout << "Variable sites present in more than " << minproportion*100 << "% of isolates will be included" << endl;
	cout << "Output will be written to " << outprefix << endl << endl;

	alignKmers(minproportion, outprefix, args, variantonly, sample);

	return 0;
}

int compareHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska compare [options] <subject split kmer files>" << endl << endl;
	cout << "Options:" << endl;
	cout << "-h\t\tPrint this help." << endl;
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << endl;
	cout << "-q <file>\tQuery split kmer file." << endl << endl;
	return 0;
}


int compareSubcommand(int argc, char *argv[]){

	vector<string> args;
	string queryfile="";

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

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
				cerr << endl << "Expecting file name after -q flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cerr << endl << "Expecting file name after -f flag" << endl << endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			cerr << endl << "Unrecognised flag " << arg << endl << endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cerr << endl << "Too few arguments" << endl;
		compareHelp();
		return 1;
	}

	if (queryfile==""){
		cerr << endl << "Query kmer file is required" << endl;
		compareHelp();
		return 1;
	}

	compareKmerFiles(queryfile, args);

	return 0;
}

int distanceHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska distance [options] <split kmer files>" << endl << endl;
	cout << "Options:" << endl;
	cout << "-c \t\tDo not print clusters files." << endl;
	cout << "-d \t\tDo not print distances file." << endl;
	cout << "-h\t\tPrint this help." << endl;
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << endl;
	cout << "-i <float>\tIdentity cutoff for defining clusters. Isolates will be \n\t\tclustered if they share at least this proportion of the \n\t\tsplit kmers in the file with fewer kmers and pass the SNP \n\t\tcutoff. [Default = 0.9]" << endl;
	cout << "-o <file>\tPrefix for output files. [Default = distances]" << endl;
	cout << "-s <int>\tSNP cutoff for defining clusters. Isolates will be clustered \n\t\tif they are separated by fewer than this number of SNPs and \n\t\tpass the identity cutoff. [Default = 20]" << endl << endl;
	return 0;
}


int distanceSubcommand(int argc, char *argv[]){

	bool distancefile=true;
	bool clusterfile=true;
	string prefix="distances";
	float minid=0.9;
	int maxsnps=20;
	vector<string> args;

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

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
				cerr << endl << "Expecting file name after -f flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-i"){
			i++;
			if (i<argc){
				try {
					minid = getfloat(argv[i]);
				}
				catch (const invalid_argument& e) {
					cerr << endl << "Expecting float between 0 and 1 after -i flag" << endl << endl;
					return 1;
				}
			}
			else {
				cerr << endl << "Expecting float between 0 and 1 after -i flag" << endl << endl;
				return 1;
			}
			if (minid < 0 || minid > 1){
				cerr << endl << "Minimum identity must be between 0 and 1" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				prefix = argv[i];
			}
			else {
				cerr << endl << "Expecting string after -o flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-s"){
			i++;
			if (i<argc){
				try {
					maxsnps = getint(argv[i]);
				}
				catch (const invalid_argument& e) {
					cerr << endl << "Expecting positive integer after -s flag" << endl << endl;
					return 1;
				}
			}
			else {
				cerr << endl << "Expecting positive integer after -s flag" << endl << endl;
				return 1;
			}
			if (maxsnps < 0){
				cerr << endl << "Maximum number of SNPs must be 0 or greater" << endl << endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			cerr << endl << "Unrecognised flag " << arg << endl << endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}
	
	if (args.size()<1){
		cerr << endl << "Too few arguments" << endl;
		distanceHelp();
		return 1;
	}

	if (distancefile==false && clusterfile==false){
		cerr << endl << "At least one output file is required. i.e. you cannot use the -c and -d flags together." << endl << endl;
		return 1;
	}
	if (prefix=="false"){
		cerr << endl << "An output prefix (-p) is required." << endl << endl;
		return 1;
	}

	skaVersion();

	cout << "Calculating pairwise distances for ";
	if (args.size()>maxNamesToPrint){
		cout << args.size() << " files";
	}
	else {
		for (auto it = args.begin(); it != args.end(); ++it){
			cout << *it << " ";
		}
	}
	cout << endl;


	cout << "Clusters will be created for isolates that are within " << maxsnps << " SNPs of one another and share at least " << minid*100 << "% of the split kmers of the isolate with fewer kmers" << endl;

	if (distancefile){
		cout << "Distances will be written to " << prefix << ".distances.tsv" << endl;
	}
	if (clusterfile){
		cout << "Clusters will be written to files with the prefix " << prefix << endl;
	}
	cout << endl;

	kmerDistance(prefix, distancefile, clusterfile, args, maxsnps, minid);

	return 0;
}

int fastaHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska fasta [options] <fasta files>" << endl << endl;
	cout << "Options:" << endl;
	cout << "-h\t\tPrint this help." << endl;
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << endl;
	cout << "-k <int>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\t[Default = 15]" << endl;
	cout << "-o <file>\tOutput file prefix. [Default = fasta]" << endl << endl;
	return 0;
}


int fastaSubcommand(int argc, char *argv[]){

	string outfile="fasta.kmers";
	long int kmersize=15;
	vector<string> args;

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			fastaHelp();
			return 0;
		}
		else if (arg == "-k"){
			i++;
			if (i<argc){
				try {
					kmersize = getint(argv[i]);
				}
				catch (const invalid_argument& e) {
					cerr << endl << "Expecting positive integer after -k flag" << endl << endl;
					return 1;
				}
			}
			else {
				cerr << endl << "Expecting positive integer after -k flag" << endl << endl;
				return 1;
			}
			if (kmersize < MinKmer || kmersize > MaxKmer || kmersize % 3 != 0){
				cerr << endl << "Kmer size must be a positive integer between " << MinKmer << " and " << MaxKmer << " and must be divisible by 3" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = string(argv[i])+".kmers";
			}
			else {
				cerr << endl << "Expecting file name after -o flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cerr << endl << "Expecting file name after -f flag" << endl << endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			cerr << endl << "Unrecognised flag " << arg << endl << endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cerr << endl << "Too few arguments" << endl;
		fastaHelp();
		return 1;
	}

	skaVersion();

	cout << "Creating split " << kmersize << "mers for ";
	if (args.size()>maxNamesToPrint){
		cout << args.size() << " files";
	}
	else {
		for (auto it = args.begin(); it != args.end(); ++it){
			cout << *it << " ";
		}
	}
	cout << endl;
	cout << "Output will be written to " << outfile << endl;

	cout << endl;

	fastaToKmers(args, outfile, kmersize);

	return 0;
}



int fastqHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska fastq [options] <fastq files>" << endl << endl;
	cout << "Options:" << endl;
	cout << "-h\t\tPrint this help." << endl;
	cout << "-c <int>\tCoverage cutoff. Kmers with coverage below this value will \n\t\tbe discarded. [Default = 4]" << endl;
	cout << "-C <int>\tFile coverage cutoff. Kmers with coverage below this value \n\t\tin any of the fastq files will be discarded. [Default = 2]" << endl;
	cout << "-k <int>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\t[Default = 15]" << endl;
	cout << "-m <float>\tMinimum allowable minor allele frequency. Kmer alleles below \n\t\tthis frequency will be discarded. [Default = 0.2]" << endl;
	cout << "-o <file>\tOutput prefix. [Default = fastq]" << endl;
	cout << "-q <int>\tQuality filter for fastq files. No kmers will be created \n\t\tfrom sequence including quality scores below this cutoff. \n\t\t[Default = 20]" << endl << endl;
	return 0;
}


int fastqSubcommand(int argc, char *argv[]){

	string outfile="fastq.kmers";
	long int covcutoff=4;
	long int filecovcutoff=2;
	long int kmersize=15;
	long int minquality=20;
	float minmaf=0.2;
	vector<string> args;

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			fastqHelp();
			return 0;
		}
		else if (arg == "-c"){
			i++;
			if (i<argc){
				try {
					covcutoff = getint(argv[i]);
				}
				catch (const invalid_argument& e) {
					cerr << endl << "Expecting positive integer after -c flag" << endl << endl;
					return 1;
				}
			}
			else {
				cerr << endl << "Expecting positive integer after -c flag" << endl << endl;
				return 1;
			}
			if (covcutoff < MinCov){
				cerr << endl << "Coverage cutoff must be a positive integer greater than or equal to " << MinCov << endl << endl;
				return 1;
			}
		}
		else if (arg == "-C"){
			i++;
			if (i<argc){
				try {
					filecovcutoff = getint(argv[i]);
				}
				catch (const invalid_argument& e) {
					cerr << endl << "Expecting positive integer after -C flag" << endl << endl;
					return 1;
				}
			}
			else {
				cerr << endl << "Expecting positive integer after -C flag" << endl << endl;
				return 1;
			}
			if (filecovcutoff <= MinCov){
				cerr << endl << "File coverage cutoff must be a positive integer greater than or equal to " << MinCov << endl << endl;
				return 1;
			}
		}
		else if (arg == "-k"){
			i++;
			if (i<argc){
				try {
					kmersize = getint(argv[i]);
				}
				catch (const invalid_argument& e) {
					cerr << endl << "Expecting positive integer after -k flag" << endl << endl;
					return 1;
				}
			}
			else {
				cerr << endl << "Expecting positive integer after -k flag" << endl << endl;
				return 1;
			}
			if (kmersize < MinKmer || kmersize > MaxKmer || kmersize % 3 != 0){
				cerr << endl << "Kmer size must be a positive integer between " << MinKmer << " and " << MaxKmer << " and must be divisible by 3" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-m"){
			i++;
			if (i<argc){
				try {
					minmaf = getfloat(argv[i]);
				}
				catch (const invalid_argument& e) {
					cerr << endl << "Expecting float between 0 and 1 after -m flag" << endl << endl;
					return 1;
				}
			}
			else {
				cerr << endl << "Expecting float between 0 and 1 after -m flag" << endl << endl;
				return 1;
			}
			if (minmaf < 0 || minmaf > 1){
				cerr << endl << "Minimum minor allele frequency must be between 0 and 1" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = string(argv[i])+".kmers";
			}
			else {
				cerr << endl << "Expecting file name after -o flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-q"){
			i++;
			if (i<argc){
				try {
					minquality = getint(argv[i]);
				}
				catch (const invalid_argument& e) {
					cerr << endl << "Expecting positive integer after -q flag" << endl << endl;
					return 1;
				}
			}
			else {
				cerr << endl << "Expecting positive integer after -q flag" << endl << endl;
				return 1;
			}
			if (minquality < MinQual || minquality > MaxQual){
				cerr << endl << "Minimum quality cutoff must be a positive integer between " << MinQual << " and " << MaxQual << endl << endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			cerr << endl << "Unrecognised flag " << arg << endl << endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cerr << endl << "Too few arguments" << endl;
		fastqHelp();
		return 1;
	}

	skaVersion();

	cout << "Creating split " << kmersize << "mers for ";
	if (args.size()>maxNamesToPrint){
		cout << args.size() << " files";
	}
	else {
		for (auto it = args.begin(); it != args.end(); ++it){
			cout << *it << " ";
		}
	}
	cout << endl;
	cout << "Filters:" << endl;
	cout << "\tCoverage cutoff = " << covcutoff << endl;
	cout << "\tPer file coverage cutoff = " << filecovcutoff << endl;
	cout << "\tMinimum minor allele frequency = " << minmaf << endl;
	cout << "\tMinimum base quality = " << minquality << endl;
	cout << "Output will be written to " << outfile << endl;

	cout << endl;

	if (fastqToKmers(args, outfile, kmersize, minquality, filecovcutoff, covcutoff, minmaf)){return 1;}

	return 0;
}


int infoHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska info [options] <split kmer files>" << endl << endl;
	cout << "Options:" << endl;
	cout << "-h\t\tPrint this help." << endl;
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << endl;
	return 0;
}


int infoSubcommand(int argc, char *argv[]){

	vector<string> args;
	bool tabulated=false;

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

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
				cerr << endl << "Expecting file name after -f flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-t"){
			tabulated=true;
		}
		else if (arg[0]=='-'){
			cerr << endl << "Unrecognised flag " << arg << endl << endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cerr << endl << "Too few arguments" << endl;
		infoHelp();
		return 1;
	}

	getKmerFileInfo(args, tabulated);

	return 0;
}

int mapHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska map [options] <split kmer files>" << endl << endl;
	cout << "Options:" << endl;
	cout << "-a <file>\tMap all bases of kmers (Default = just map middle base)." << endl;
	cout << "-h\t\tPrint this help." << endl;
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << endl;
	cout << "-k <int>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\tMust be the same value used to create the kmer files. \n\t\t[Default = 15]" << endl;
	cout << "-i\t\tInclude reference sequence in alignment." << endl;
	cout << "-m\t\tMap bases to repeats rather than making them N." << endl;
	cout << "-o <file>\tOutput file prefix. [Default = mappedkmers]" << endl;
	cout << "-r <file>\tReference fasta file name. [Required]" << endl;
	cout << "-s <file>\tFile of sample names to include in the alignment." << endl;
	cout << "-v\tOutput variant only alignment. [Default = all sites]" << endl << endl;
	return 0;
}


int mapSubcommand(int argc, char *argv[]){

	string outprefix = "mappedkmers";
	string reference = "";
	long int kmersize=15;
	bool includeref=false;
	bool maprepeats=false;
	bool fillall=false;
	bool variantonly=false;
	vector <string> sample;
	vector<string> args;

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			mapHelp();
			return 0;
		}
		else if (arg == "-k"){
			i++;
			if (i<argc){
				try {
					kmersize = getint(argv[i]);
				}
				catch (const invalid_argument& e) {
					cerr << endl << "Expecting positive integer after -k flag" << endl << endl;
					return 1;
				}
			}
			else {
				cerr << endl << "Expecting positive integer after -k flag" << endl << endl;
				return 1;
			}
			if (kmersize < MinKmer || kmersize > MaxKmer || kmersize % 3 != 0){
				cerr << endl << "Kmer size must be a positive integer between " << MinKmer << " and " << MaxKmer << " and must be divisible by 3" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-r"){
			i++;
			if (i<argc){
				reference = argv[i];
			}
			else {
				cerr << endl << "Expecting file name after -r flag" << endl << endl;
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
				outprefix = string(argv[i]);
			}
			else {
				cerr << endl << "Expecting file name after -o flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cerr << endl << "Expecting file name after -f flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-s"){
			i++;
			if (i<argc){
				fileToVector(argv[i], sample);
			}
			else {
				cerr << endl << "Expecting file name after -s flag" << endl << endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			cerr << endl << "Unrecognised flag " << arg << endl << endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cerr << endl << "Too few arguments" << endl;
		mapHelp();
		return 1;
	}

	if (variantonly){
		outprefix+="_variants";
	}

	outprefix+=".aln";

	if (reference==""){
		cerr << endl << "Reference file is required" << endl;
		mapHelp();
		return 1;
	}

	skaVersion();

	cout << "Aligning ";
	if (args.size()>maxNamesToPrint){
		cout << args.size() << " files";
	}
	else {
		for (auto it = args.begin(); it != args.end(); ++it){
			cout << *it << " ";
		}
	}
	cout << endl;
	cout << "Using split kmer size of " << kmersize << endl;
	cout << "Using " << reference << " as reference" << endl;
	if (includeref){
		cout << "Reference will be included in the alignment" << endl;
	}
	else{
		cout << "Reference will not be included in the alignment" << endl;
	}
	if (maprepeats){
		cout << "Repeats in the reference will be mapped rather than left as gaps" << endl;
	}
	else{
		cout << "Repeats in the reference will be output as gaps" << endl;
	}
	if (fillall){
		cout << "All bases in kmers will be mapped as lower case" << endl;
	}
	cout << "Output will be written to " << outprefix << endl << endl;

	alignKmersToReference(reference, outprefix, args, kmersize, includeref, maprepeats, fillall, variantonly, sample);

	return 0;
}


int mergeHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska merge [options] <subject split kmer files>" << endl << endl;
	cout << "Options:" << endl;
	cout << "-h\t\tPrint this help." << endl;
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << endl;
	cout << "-o <file>\tOutput file prefix. [Default = merged]" << endl;
	cout << "-s <file>\tFile of sample names to include in the merged file. Can \n\t\tbe used to extract samples from a merge" << endl << endl;
	return 0;
}


int mergeSubcommand(int argc, char *argv[]){

	vector<string> args;
	vector <string> sample;
	string outfile="merged";

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			mergeHelp();
			return 0;
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = string(argv[i])+"";
			}
			else {
				cerr << endl << "Expecting file name after -o flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cerr << endl << "Expecting file name after -f flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-s"){
			i++;
			if (i<argc){
				fileToVector(argv[i], sample);
			}
			else {
				cerr << endl << "Expecting file name after -s flag" << endl << endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			cerr << endl << "Unrecognised flag " << arg << endl << endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cerr << endl << "Too few arguments" << endl;
		mergeHelp();
		return 1;
	}

	if (outfile==""){
		cerr << endl << "Output file name is required" << endl;
		mergeHelp();
		return 1;
	}

	skaVersion();

	cout << "Merging ";
	if (args.size()>maxNamesToPrint){
		cout << args.size() << " files";
	}
	else {
		for (auto it = args.begin(); it != args.end(); ++it){
			cout << *it << " ";
		}
	}
	cout << endl;
	cout << "Output will be written to " << outfile << endl << endl;

	mergeKmerFiles(outfile, args, sample);

	return 0;
}


int summaryHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska summary [options] <split kmer files>" << endl << endl;
	cout << "Options:" << endl;
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << endl;
	cout << "-h\t\tPrint this help." << endl << endl;
	return 0;
}


int summarySubcommand(int argc, char *argv[]){

	vector<string> args;

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

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
				cerr << endl << "Expecting file name after -f flag" << endl << endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			cerr << endl << "Unrecognised flag " << arg << endl << endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cerr << endl << "Too few arguments" << endl;
		compareHelp();
		return 1;
	}

	summariseKmerFiles(args);

	return 0;
}


int typeHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska type [options] <subject split kmer files>" << endl << endl;
	cout << "Options:" << endl;
	cout << "-h\t\tPrint this help." << endl;
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << endl;
	cout << "-p <file>\ttab file containing profile information." << endl;
	cout << "-q <file>\tQuery split kmer file. This can be a single kmer file or a \n\t\tkmerge." << endl << endl;
	return 0;
}


int typeSubcommand(int argc, char *argv[]){

	vector<string> args;
	string queryfile="";
	string profilefile="";

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

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
				cerr << endl << "Expecting file name after -q flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cerr << endl << "Expecting file name after -f flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-p"){
			i++;
			if (i<argc){
				profilefile = argv[i];
			}
			else {
				cerr << endl << "Expecting file name after -p flag" << endl << endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			cerr << endl << "Unrecognised flag " << arg << endl << endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cerr << endl << "Too few arguments" << endl;
		typeHelp();
		return 1;
	}

	if (queryfile==""){
		cerr << endl << "Query kmer file is required" << endl;
		typeHelp();
		return 1;
	}

	typeKmerFile(queryfile, profilefile, args);

	return 0;
}


int uniqueHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska unique [options]" << endl << endl;
	cout << "Options:" << endl;
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << endl;
	cout << "-h\t\tPrint this help." << endl;
	cout << "-i <file>\tFile of ingroup sample names. Unique kmers found \n\t\tin these files will be retained." << endl;
	cout << "-o <file>\tOutput file prefix. [Default = unique]" << endl;
	cout << "-p <float>\tMinimum proportion of ingroup isolates required to possess a \n\t\tsplit kmer for that kmer to be retained. [Default = 0.9]" << endl << endl;
					   //123456789012345678901234567890123456789012345678901234567890
	return 0;
}


int uniqueSubcommand(int argc, char *argv[]){

	string outfile="unique.kmers";
	vector<string> ingroup;
	float minproportion=0.9;
	vector<string> args;

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			uniqueHelp();
			return 0;
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = string(argv[i])+".kmers";
			}
			else {
				cerr << endl << "Expecting file name after -o flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-i"){
			i++;
			if (i<argc){
				fileToVector(argv[i], ingroup);
			}
			else {
				cerr << endl << "Expecting file name after -i flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cerr << endl << "Expecting file name after -f flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-p"){
			i++;
			if (i<argc){
				try {
					minproportion = getfloat(argv[i]);
				}
				catch (const invalid_argument& e) {
					cerr << endl << "Expecting float between 0 and 1 after -p flag" << endl << endl;
					return 1;
				}
			}
			else {
				cerr << endl << "Expecting float between 0 and 1 after -p flag" << endl << endl;
				return 1;
			}
			if (minproportion < 0 || minproportion > 1){
				cerr << endl << "Minimum proportion must be between 0 and 1" << endl << endl;
				return 1;
			}
		}
		else {
			args.push_back(arg);
		}
	}

	if (ingroup.size()<1){
		cerr << endl << "At least one ingroup kmer file must be specified" << endl;
		uniqueHelp();
		return 1;
	}

	if (args.size()<1){
		cerr << endl << "At least one kmer file must be specified" << endl;
		uniqueHelp();
		return 1;
	}

	skaVersion();

	cout << "Finding unique kmers in at least " << minproportion*100 << "% of ";
	for (auto it = ingroup.begin(); it != ingroup.end(); ++it){
		cout << *it << " ";
	}
	/*cout << endl << "That are not also found in ";
	for (auto it = outgroup.begin(); it != args.end(); ++it){
		cout << *it << " ";
	}*/
	cout << endl;

	cout << "Output will be written to " << outfile << endl;

	cout << endl;

	uniqueKmers(ingroup, args, minproportion, outfile);

	return 0;
}

int weedHelp(void){
	cout << endl << "Usage:" << endl;
	cout << "ska weed [options] <split kmer files>" << endl << endl;
	cout << "Options:" << endl;
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line." << endl;
	cout << "-h\t\tPrint this help." << endl;
	cout << "-i <file>\tName of kmer file containing kmers to be weeded. [Required]" << endl;
	return 0;
}


int weedSubcommand(int argc, char *argv[]){

	string infile="";
	string kmerfile="";
	vector<string> args;

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			weedHelp();
			return 0;
		}
		else if (arg == "-i"){
			i++;
			if (i<argc){
				infile = argv[i];
			}
			else {
				cerr << endl << "Expecting file name after -i flag" << endl << endl;
				return 1;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cerr << endl << "Expecting file name after -f flag" << endl << endl;
				return 1;
			}
		}
		else if (arg[0]=='-'){
			cerr << endl << "Unrecognised flag " << arg << endl << endl;
				return 1;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cerr << endl << "Too few arguments" << endl;
		weedHelp();
		return 1;
	}

	if (infile==""){
		cerr << "Input split kmer file is required" << endl;
		weedHelp();
		return 1;
	}

	skaVersion();

	cout << "Weeding kmers in " << infile << endl;
	cout << "From ";
	if (args.size()>maxNamesToPrint){
		cout << args.size() << " files";
	}
	else {
		for (auto it = args.begin(); it != args.end(); ++it){
			cout << *it << " ";
		}
	}

	cout << endl;

	weedKmers(args, infile);

	return 0;
}


int main(int argc, char *argv[])
{
	
	if (argc<2){
		skaHelp();
		return 0;
	}

	string subcommand=(argv[1]);

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
			cerr << endl << "Unrecognised flag " << subcommand << endl << endl;
			return 1;
		}
	else {
		cerr << endl << "Unrecognised subcommand" << endl << endl;
		skaHelp();
		return 1;
	}
	
	return 0;

}
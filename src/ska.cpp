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
#include "kmers.hpp"

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
	cout << "Usage:\n";
	cout << "ska <subcommand>\n\n";
	cout << "Subcommands:\n";
	cout << "align\t\tReference-free alignment from a set of split kmer files\n";
	cout << "alleles\t\tCreate a merged split kmer file for all sequenes in one or \n\t\tmore multifasta files\n";
	cout << "compare\t\tPrint comparison statistics for a query split kmer file \n\t\tagainst a set of subject split kmer files\n";
	cout << "distance\tPairwise distance calculation and clustering from split kmer \n\t\tfiles\n";
	cout << "fasta\t\tCreate a split kmer file from fasta file(s)\n";
	cout << "fastq\t\tCreate a split kmer file from fastq file(s)\n";
	cout << "info\t\tPrint some information about one or more kmer or kmerge \n\t\tfiles\n";
				   //123456789012345678901234567890123456789012345678901234567890
	cout << "map\t\tAlign split kmer file(s) against a reference fasta file\n";
	cout << "summary\t\tPrint split kmer file summary statistics\n";
	cout << "type\t\tType split kmer files using a set of allele files\n";
	cout << "unique\t\tOutput kmers unique to a set of split kmer files\n";
	cout << "version\t\tPrint the version and citation for ska\n";
	cout << "weed\t\tWeed kmers from a split kmer file\n\n";

	               
	return 0;
}

int skaVersion(void){
	cout << "\nSKA: Split Kmer Analysis\n";
	cout << "Version: " << versionNumber << "\n";
	cout << "Citation: " << citation << "\n\n";
	return 0;
}
int allelesHelp(void){
	cout << "\nUsage:\n";
	cout << "ska fasta [options] <fasta files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help.\n";
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line.\n";
	cout << "-k <int>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\t[Default = 15]\n\n";
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
					cout << "\nExpecting positive integer after -k flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting positive integer after -k flag\n\n";
				return 0;
			}
			if (kmersize < MinKmer || kmersize > MaxKmer || kmersize % 3 != 0){
				cout << "\nKmer size must be a positive integer between " << MinKmer << " and " << MaxKmer << " and must be divisible by 3\n\n";
				return 0;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cout << "\nExpecting file name after -f flag\n\n";
				return 0;
			}
		}
		else if (arg[0]=='-'){
			cout << "\nUnrecognised flag " << arg << "\n\n";
				return 0;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cout << "\nToo few arguments\n";
		allelesHelp();
		return 0;
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

	cout << "\n";

	allelesToKmers(args, kmersize);

	return 0;
}



int alignHelp(void){
	cout << "\nUsage:\n";
	cout << "ska align [options] <split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help.\n";
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line.\n";
	cout << "-o <file>\tPrefix for output files. [Default = reference_free]\n";
	cout << "-p <float>\tMinimum proportion of isolates required to possess a split \n\t\tkmer for that kmer to be included in the alignment. \n\t\t[Default = 0.9]\n";
	cout << "-s <file>\tFile of sample names to include in the alignment.\n";
	cout << "-v\tOutput variant only alignment. [Default = all sites]\n\n";
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
				cout << "\nExpecting string after -o flag\n\n";
				return 0;
			}
		}
		else if (arg == "-p"){
			i++;
			if (i<argc){
				try {
					minproportion = getfloat(argv[i]);
				}
				catch (const invalid_argument& e) {
					cout << "\nExpecting float between 0 and 1 after -p flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting float between 0 and 1 after -p flag\n\n";
				return 0;
			}
			if (minproportion < 0 || minproportion > 1){
				cout << "\nMinimum proportion must be between 0 and 1\n\n";
				return 0;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cout << "\nExpecting file name after -f flag\n\n";
				return 0;
			}
		}
		else if (arg == "-s"){
			i++;
			if (i<argc){
				fileToVector(argv[i], sample);
			}
			else {
				cout << "\nExpecting file name after -s flag\n\n";
				return 0;
			}
		}
		else if (arg == "-v"){
			variantonly = true;
		}
		else if (arg[0]=='-'){
			cout << "\nUnrecognised flag " << arg << "\n\n";
				return 0;
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
		cout << "\nToo few arguments\n";
		alignHelp();
		return 0;
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
	cout << "\n";
	cout << "Variable sites present in more than " << minproportion*100 << "% of isolates will be included\n";
	cout << "Output will be written to " << outprefix << "\n\n";

	alignKmers(minproportion, outprefix, args, variantonly, sample);

	return 0;
}

int compareHelp(void){
	cout << "\nUsage:\n";
	cout << "ska compare [options] <subject split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help.\n";
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line.\n";
	cout << "-q <file>\tQuery split kmer file.\n\n";
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
				cout << "\nExpecting file name after -q flag\n\n";
				return 0;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cout << "\nExpecting file name after -f flag\n\n";
				return 0;
			}
		}
		else if (arg[0]=='-'){
			cout << "\nUnrecognised flag " << arg << "\n\n";
				return 0;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cout << "\nToo few arguments\n";
		compareHelp();
		return 0;
	}

	if (queryfile==""){
		cout << "\nQuery kmer file is required\n";
		compareHelp();
		return 0;
	}

	compareKmerFiles(queryfile, args);

	return 0;
}

int distanceHelp(void){
	cout << "\nUsage:\n";
	cout << "ska distance [options] <split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-c \t\tDo not print clusters files.\n";
	cout << "-d \t\tDo not print distances file.\n";
	cout << "-h\t\tPrint this help.\n";
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line.\n";
	cout << "-i <float>\tIdentity cutoff for defining clusters. Isolates will be \n\t\tclustered if they share at least this proportion of the \n\t\tsplit kmers in the file with fewer kmers and pass the SNP \n\t\tcutoff. [Default = 0.9]\n";
	cout << "-o <file>\tPrefix for output files. [Default = distances]\n";
	cout << "-s <int>\tSNP cutoff for defining clusters. Isolates will be clustered \n\t\tif they are separated by fewer than this number of SNPs and \n\t\tpass the identity cutoff. [Default = 20]\n\n";
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
				cout << "\nExpecting file name after -f flag\n\n";
				return 0;
			}
		}
		else if (arg == "-i"){
			i++;
			if (i<argc){
				try {
					minid = getfloat(argv[i]);
				}
				catch (const invalid_argument& e) {
					cout << "\nExpecting float between 0 and 1 after -i flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting float between 0 and 1 after -i flag\n\n";
				return 0;
			}
			if (minid < 0 || minid > 1){
				cout << "\nMinimum identity must be between 0 and 1\n\n";
				return 0;
			}
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				prefix = argv[i];
			}
			else {
				cout << "\nExpecting string after -o flag\n\n";
				return 0;
			}
		}
		else if (arg == "-s"){
			i++;
			if (i<argc){
				try {
					maxsnps = getint(argv[i]);
				}
				catch (const invalid_argument& e) {
					cout << "\nExpecting positive integer after -s flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting positive integer after -s flag\n\n";
				return 0;
			}
			if (maxsnps < 0){
				cout << "\nMaximum number of SNPs must be 0 or greater\n\n";
				return 0;
			}
		}
		else if (arg[0]=='-'){
			cout << "\nUnrecognised flag " << arg << "\n\n";
				return 0;
		}
		else {
			args.push_back(arg);
		}
	}
	
	if (args.size()<1){
		cout << "\nToo few arguments\n";
		distanceHelp();
		return 0;
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
	cout << "\n";


	cout << "Clusters will be created for isolates that are within " << maxsnps << " SNPs of one another and share at least " << minid*100 << "% of the split kmers of the isolate with fewer kmers\n";

	if (distancefile){
		cout << "Distances will be written to " << prefix << ".distances.tsv\n";
	}
	if (clusterfile){
		cout << "Clusters will be written to files with the prefix " << prefix << "\n";
	}
	if (distancefile==false && clusterfile==false){
		cout << "\nAt least one output file is required. i.e. you cannot use the -c and -d flags together.\n\n";
		return 0;
	}
	if (prefix=="false"){
		cout << "\nAn output prefix (-p) is required.\n\n";
		return 0;
	}
	cout << "\n";

	kmerDistance(prefix, distancefile, clusterfile, args, maxsnps, minid);

	return 0;
}

int fastaHelp(void){
	cout << "\nUsage:\n";
	cout << "ska fasta [options] <fasta files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help.\n";
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line.\n";
	cout << "-k <int>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\t[Default = 15]\n";
	cout << "-o <file>\tOutput file prefix. [Default = fasta]\n\n";
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
					cout << "\nExpecting positive integer after -k flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting positive integer after -k flag\n\n";
				return 0;
			}
			if (kmersize < MinKmer || kmersize > MaxKmer || kmersize % 3 != 0){
				cout << "\nKmer size must be a positive integer between " << MinKmer << " and " << MaxKmer << " and must be divisible by 3\n\n";
				return 0;
			}
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = string(argv[i])+".kmers";
			}
			else {
				cout << "\nExpecting file name after -o flag\n\n";
				return 0;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cout << "\nExpecting file name after -f flag\n\n";
				return 0;
			}
		}
		else if (arg[0]=='-'){
			cout << "\nUnrecognised flag " << arg << "\n\n";
				return 0;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cout << "\nToo few arguments\n";
		fastaHelp();
		return 0;
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

	cout << "\n";

	fastaToKmers(args, outfile, kmersize);

	return 0;
}



int fastqHelp(void){
	cout << "\nUsage:\n";
	cout << "ska fastq [options] <fastq files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help.\n";
	cout << "-c <int>\tCoverage cutoff. Kmers with coverage below this value will \n\t\tbe discarded. [Default = 4]\n";
	cout << "-C <int>\tFile coverage cutoff. Kmers with coverage below this value \n\t\tin any of the fastq files will be discarded. [Default = 2]\n";
	cout << "-k <int>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\t[Default = 15]\n";
	cout << "-m <float>\tMinimum allowable minor allele frequency. Kmer alleles below \n\t\tthis frequency will be discarded. [Default = 0.2]\n";
	cout << "-o <file>\tOutput prefix. [Default = fastq]\n";
	cout << "-q <int>\tQuality filter for fastq files. No kmers will be created \n\t\tfrom sequence including quality scores below this cutoff. \n\t\t[Default = 20]\n\n";
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
					cout << "\nExpecting positive integer after -c flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting positive integer after -c flag\n\n";
				return 0;
			}
			if (covcutoff < MinCov){
				cout << "\nCoverage cutoff must be a positive integer greater than or equal to " << MinCov << "\n\n";
				return 0;
			}
		}
		else if (arg == "-C"){
			i++;
			if (i<argc){
				try {
					filecovcutoff = getint(argv[i]);
				}
				catch (const invalid_argument& e) {
					cout << "\nExpecting positive integer after -C flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting positive integer after -C flag\n\n";
				return 0;
			}
			if (filecovcutoff <= MinCov){
				cout << "\nFile coverage cutoff must be a positive integer greater than or equal to " << MinCov << "\n\n";
				return 0;
			}
		}
		else if (arg == "-k"){
			i++;
			if (i<argc){
				try {
					kmersize = getint(argv[i]);
				}
				catch (const invalid_argument& e) {
					cout << "\nExpecting positive integer after -k flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting positive integer after -k flag\n\n";
				return 0;
			}
			if (kmersize < MinKmer || kmersize > MaxKmer || kmersize % 3 != 0){
				cout << "\nKmer size must be a positive integer between " << MinKmer << " and " << MaxKmer << " and must be divisible by 3\n\n";
				return 0;
			}
		}
		else if (arg == "-m"){
			i++;
			if (i<argc){
				try {
					minmaf = getfloat(argv[i]);
				}
				catch (const invalid_argument& e) {
					cout << "\nExpecting float between 0 and 1 after -m flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting float between 0 and 1 after -m flag\n\n";
				return 0;
			}
			if (minmaf < 0 || minmaf > 1){
				cout << "\nMinimum minor allele frequency must be between 0 and 1\n\n";
				return 0;
			}
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = string(argv[i])+".kmers";
			}
			else {
				cout << "\nExpecting file name after -o flag\n\n";
				return 0;
			}
		}
		else if (arg == "-q"){
			i++;
			if (i<argc){
				try {
					minquality = getint(argv[i]);
				}
				catch (const invalid_argument& e) {
					cout << "\nExpecting positive integer after -q flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting positive integer after -q flag\n\n";
				return 0;
			}
			if (minquality < MinQual || minquality > MaxQual){
				cout << "\nMinimum quality cutoff must be a positive integer between " << MinQual << " and " << MaxQual << "\n\n";
				return 0;
			}
		}
		else if (arg[0]=='-'){
			cout << "\nUnrecognised flag " << arg << "\n\n";
				return 0;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cout << "\nToo few arguments\n";
		fastqHelp();
		return 0;
	}
	else if (args.size()>2){
		cout << "\nA maximum of two fastq files can be used\n";
		fastqHelp();
		return 0;
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
	cout << "\n";
	cout << "Filters:\n";
	cout << "\tCoverage cutoff = " << covcutoff << endl;
	cout << "\tPer file coverage cutoff = " << filecovcutoff << endl;
	cout << "\tMinimum minor allele frequency = " << minmaf << endl;
	cout << "\tMinimum base quality = " << minquality << endl;
	cout << "Output will be written to " << outfile << endl;

	cout << "\n";

	fastqToKmers(args, outfile, kmersize, minquality, filecovcutoff, covcutoff, minmaf);

	return 0;
}


int infoHelp(void){
	cout << "\nUsage:\n";
	cout << "ska info [options] <split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help.\n";
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line.\n";
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
				cout << "\nExpecting file name after -f flag\n\n";
				return 0;
			}
		}
		else if (arg == "-t"){
			tabulated=true;
		}
		else if (arg[0]=='-'){
			cout << "\nUnrecognised flag " << arg << "\n\n";
				return 0;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cout << "\nToo few arguments\n";
		infoHelp();
		return 0;
	}

	getKmerFileInfo(args, tabulated);

	return 0;
}

int mapHelp(void){
	cout << "\nUsage:\n";
	cout << "ska map [options] <split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-a <file>\tMap all bases of kmers (Default = just map middle base).\n";
	cout << "-h\t\tPrint this help.\n";
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line.\n";
	cout << "-k <int>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\tMust be the same value used to create the kmer files. \n\t\t[Default = 15]\n";
	cout << "-i\t\tInclude reference sequence in alignment.\n";
	cout << "-m\t\tMap bases to repeats rather than making them N.\n";
	cout << "-o <file>\tOutput file prefix. [Default = mappedkmers]\n";
	cout << "-r <file>\tReference fasta file name. [Required]\n";
	cout << "-s <file>\tFile of sample names to include in the alignment.\n";
	cout << "-v\tOutput variant only alignment. [Default = all sites]\n\n";
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
					cout << "\nExpecting positive integer after -k flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting positive integer after -k flag\n\n";
				return 0;
			}
			if (kmersize < MinKmer || kmersize > MaxKmer || kmersize % 3 != 0){
				cout << "\nKmer size must be a positive integer between " << MinKmer << " and " << MaxKmer << " and must be divisible by 3\n\n";
				return 0;
			}
		}
		else if (arg == "-r"){
			i++;
			if (i<argc){
				reference = argv[i];
			}
			else {
				cout << "\nExpecting file name after -r flag\n\n";
				return 0;
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
				cout << "\nExpecting file name after -o flag\n\n";
				return 0;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cout << "\nExpecting file name after -f flag\n\n";
				return 0;
			}
		}
		else if (arg == "-s"){
			i++;
			if (i<argc){
				fileToVector(argv[i], sample);
			}
			else {
				cout << "\nExpecting file name after -s flag\n\n";
				return 0;
			}
		}
		else if (arg[0]=='-'){
			cout << "\nUnrecognised flag " << arg << "\n\n";
				return 0;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cout << "\nToo few arguments\n";
		mapHelp();
		return 0;
	}

	if (variantonly){
		outprefix+="_variants";
	}

	outprefix+=".aln";

	if (reference==""){
		cout << "\nReference file is required\n";
		mapHelp();
		return 0;
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
	cout << "\n";
	cout << "Using split kmer size of " << kmersize << "\n";
	cout << "Using " << reference << " as reference\n";
	if (includeref){
		cout << "Reference will be included in the alignment\n";
	}
	else{
		cout << "Reference will not be included in the alignment\n";
	}
	if (maprepeats){
		cout << "Repeats in the reference will be mapped rather than left as gaps\n";
	}
	else{
		cout << "Repeats in the reference will be output as gaps\n";
	}
	if (fillall){
		cout << "All bases in kmers will be mapped as lower case\n";
	}
	cout << "Output will be written to " << outprefix << "\n\n";

	alignKmersToReference(reference, outprefix, args, kmersize, includeref, maprepeats, fillall, variantonly, sample);

	return 0;
}


int mergeHelp(void){
	cout << "\nUsage:\n";
	cout << "ska merge [options] <subject split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help.\n";
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line.\n";
	cout << "-o <file>\tOutput file prefix. [Default = merged]\n";
	cout << "-s <file>\tFile of sample names to include in the merged file. Can \n\t\tbe used to extract samples from a merge\n\n";
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
				cout << "\nExpecting file name after -o flag\n\n";
				return 0;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cout << "\nExpecting file name after -f flag\n\n";
				return 0;
			}
		}
		else if (arg == "-s"){
			i++;
			if (i<argc){
				fileToVector(argv[i], sample);
			}
			else {
				cout << "\nExpecting file name after -s flag\n\n";
				return 0;
			}
		}
		else if (arg[0]=='-'){
			cout << "\nUnrecognised flag " << arg << "\n\n";
				return 0;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cout << "\nToo few arguments\n";
		mergeHelp();
		return 0;
	}

	if (outfile==""){
		cout << "\nOutput file name is required\n";
		mergeHelp();
		return 0;
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
	cout << "\n";
	cout << "Output will be written to " << outfile << "\n\n";

	mergeKmerFiles(outfile, args, sample);

	return 0;
}


int summaryHelp(void){
	cout << "\nUsage:\n";
	cout << "ska summary [options] <split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line.\n";
	cout << "-h\t\tPrint this help.\n\n";
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
				cout << "\nExpecting file name after -f flag\n\n";
				return 0;
			}
		}
		else if (arg[0]=='-'){
			cout << "\nUnrecognised flag " << arg << "\n\n";
				return 0;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cout << "\nToo few arguments\n";
		compareHelp();
		return 0;
	}

	summariseKmerFiles(args);

	return 0;
}


int typeHelp(void){
	cout << "\nUsage:\n";
	cout << "ska type [options] <subject split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help.\n";
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line.\n";
	cout << "-p <file>\ttab file containing profile information.\n";
	cout << "-q <file>\tQuery split kmer file. This can be a single kmer file or a \n\t\tkmerge.\n\n";
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
				cout << "\nExpecting file name after -q flag\n\n";
				return 0;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cout << "\nExpecting file name after -f flag\n\n";
				return 0;
			}
		}
		else if (arg == "-p"){
			i++;
			if (i<argc){
				profilefile = argv[i];
			}
			else {
				cout << "\nExpecting file name after -p flag\n\n";
				return 0;
			}
		}
		else if (arg[0]=='-'){
			cout << "\nUnrecognised flag " << arg << "\n\n";
				return 0;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cout << "\nToo few arguments\n";
		typeHelp();
		return 0;
	}

	if (queryfile==""){
		cout << "\nQuery kmer file is required\n";
		typeHelp();
		return 0;
	}

	typeKmerFile(queryfile, profilefile, args);

	return 0;
}


int uniqueHelp(void){
	cout << "\nUsage:\n";
	cout << "ska unique [options]\n\n";
	cout << "Options:\n";
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line.\n";
	cout << "-h\t\tPrint this help.\n";
	cout << "-i <file>\tFile of ingroup sample names. Unique kmers found \n\t\tin these files will be retained.\n";
	cout << "-o <file>\tOutput file prefix. [Default = unique]\n";
	cout << "-p <float>\tMinimum proportion of ingroup isolates required to possess a \n\t\tsplit kmer for that kmer to be retained. [Default = 0.9]\n\n";
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
				cout << "\nExpecting file name after -o flag\n\n";
				return 0;
			}
		}
		else if (arg == "-i"){
			i++;
			if (i<argc){
				fileToVector(argv[i], ingroup);
			}
			else {
				cout << "\nExpecting file name after -i flag\n\n";
				return 0;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cout << "\nExpecting file name after -f flag\n\n";
				return 0;
			}
		}
		else if (arg == "-p"){
			i++;
			if (i<argc){
				try {
					minproportion = getfloat(argv[i]);
				}
				catch (const invalid_argument& e) {
					cout << "\nExpecting float between 0 and 1 after -p flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting float between 0 and 1 after -p flag\n\n";
				return 0;
			}
			if (minproportion < 0 || minproportion > 1){
				cout << "\nMinimum proportion must be between 0 and 1\n\n";
				return 0;
			}
		}
		else {
			args.push_back(arg);
		}
	}

	if (ingroup.size()<1){
		cout << "\nAt least one ingroup kmer file must be specified\n";
		uniqueHelp();
		return 0;
	}

	if (args.size()<1){
		cout << "\nAt least one kmer file must be specified\n";
		uniqueHelp();
		return 0;
	}

	skaVersion();

	cout << "Finding unique kmers in at least " << minproportion*100 << "% of ";
	for (auto it = ingroup.begin(); it != ingroup.end(); ++it){
		cout << *it << " ";
	}
	/*cout << "\nThat are not also found in ";
	for (auto it = outgroup.begin(); it != args.end(); ++it){
		cout << *it << " ";
	}*/
	cout << "\n";

	cout << "Output will be written to " << outfile << "\n";

	cout << "\n";

	uniqueKmers(ingroup, args, minproportion, outfile);

	return 0;
}

int weedHelp(void){
	cout << "\nUsage:\n";
	cout << "ska weed [options] <split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-f <file>\tFile of split kmer file names. These will be added to or \n\t\tused as an alternative input to the list provided on the \n\t\tcommand line.\n";
	cout << "-h\t\tPrint this help.\n";
	cout << "-i <file>\tName of kmer file containing kmers to be weeded. [Required]\n";
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
				cout << "\nExpecting file name after -i flag\n\n";
				return 0;
			}
		}
		else if (arg == "-f"){
			i++;
			if (i<argc){
				fileToVector(argv[i], args);
			}
			else {
				cout << "\nExpecting file name after -f flag\n\n";
				return 0;
			}
		}
		else if (arg[0]=='-'){
			cout << "\nUnrecognised flag " << arg << "\n\n";
				return 0;
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()<1){
		cout << "\nToo few arguments\n";
		weedHelp();
		return 0;
	}

	if (infile==""){
		cout << "Input split kmer file is required\n";
		weedHelp();
		return 0;
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
			allelesSubcommand(argc, argv);
		}
	else if (subcommand == "align"){
			alignSubcommand(argc, argv);
		}
	else if (subcommand == "compare"){
			compareSubcommand(argc, argv);
		}
	else if (subcommand == "distance"){
			distanceSubcommand(argc, argv);
		}
	else if (subcommand == "fasta"){
			fastaSubcommand(argc, argv);
		}
	else if (subcommand == "fastq"){
			fastqSubcommand(argc, argv);
		}
	else if (subcommand == "info"){
			infoSubcommand(argc, argv);
		}
	else if (subcommand == "map"){
			mapSubcommand(argc, argv);
		}
	else if (subcommand == "merge"){
			mergeSubcommand(argc, argv);
		}
	else if (subcommand == "summary"){
			summarySubcommand(argc, argv);
		}
	else if (subcommand == "type"){
			typeSubcommand(argc, argv);
		}
	else if (subcommand == "unique"){
			uniqueSubcommand(argc, argv);
		}
	else if (subcommand == "weed"){
			weedSubcommand(argc, argv);
		}
	else if (subcommand == "-v" || subcommand == "--version" || subcommand == "version"){
			skaVersion();
		}
	else if (subcommand == "-h" || subcommand == "--help"){
			skaHelp();
		}
	else if (subcommand[0]=='-'){
			cout << "\nUnrecognised flag " << subcommand << "\n\n";
			return 0;
		}
	else {
		cout << "\nUnrecognised subcommand\n\n";
		skaHelp();
	}
	
	return 0;

}
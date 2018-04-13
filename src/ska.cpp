#include <iostream> //std::cout
#include <string> //std::string
#include <vector> //std::string
#include <stdlib.h> //std::strtol
#include <stdexcept>      // std::invalid_argument
#include "sk_align.hpp"
#include "sk_compare.hpp"
#include "sk_fasta.hpp"
#include "sk_fastq.hpp"
#include "sk_refalign.hpp"
#include "sk_summary.hpp"
#include "sk_weed.hpp"

using namespace std;

long int getint(char * arg);
float getfloat(char * arg);
int skaHelp(void);
int alignHelp(void);
int alignSubcommand(int argc, char *argv[]);
int compareHelp(void);
int compareSubcommand(int argc, char *argv[]);
int fastaHelp(void);
int fastaSubcommand(int argc, char *argv[]);
int fastqHelp(void);
int fastqSubcommand(int argc, char *argv[]);
int mapHelp(void);
int mapSubcommand(int argc, char *argv[]);
int summaryHelp(void);
int summarySubcommand(int argc, char *argv[]);
int weedHelp(void);
int weedSubcommand(int argc, char *argv[]);

int MinKmer=9;
int MaxKmer=30;
int MinQual=0;
int MaxQual=60;
int MinCov=0;

string versionNumber = "0.1";

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
	cout << "\nSKA: Split Kmer Analysis\n";
	cout << "\nUsage:\n";
	cout << "ska <subcommand>\n\n";
	cout << "Subcommands:\n";
	cout << "align\t\tReference-free alignment from split kmer files\n";
	cout << "compare\t\tCompare two split kmer files\n";
	cout << "fasta\t\tCreate split kmer file from fasta file(s)\n";
	cout << "fastq\t\tCreate split kmer file from fastq file(s)\n";
	cout << "map\t\tAlign split kmer file(s) against a reference fasta file\n";
	cout << "summary\t\tPrint split kmer file summary statistics\n";
	cout << "version\t\tPrint the version and citation for ska\n";
	cout << "weed\t\tWeed kmers from a split kmer file\n\n";
	return 0;
}

int skaVersion(void){
	cout << "\nSKA: Split Kmer Analysis\n";
	cout << "Version: " << versionNumber << "\n";
	cout << "Citation: TBA\n\n";
	return 0;
}

int alignHelp(void){
	cout << "\nUsage:\n";
	cout << "ska align [options] <split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help\n";
	cout << "-o\t\tOutput file name [Default = reference_free.aln]\n";
	cout << "-p <float>\tMinimum proportion of isolates required to possess a split \n\t\tkmer for that kmer to be included in the alignment. \n\t\t[Default = 0.1]\n\n";
	return 0;
}


int alignSubcommand(int argc, char *argv[]){

	string outfile="reference_free.aln";
	float minproportion=0.1;
	vector<string> args;

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			alignHelp();
			return 0;
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = argv[i];
			}
			else {
				cout << "\nExpecting file name after -o flag\n\n";
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
	
	if (args.size()<2){
		cout << "\nToo few arguments\n";
		alignHelp();
		return 0;
	}

	cout << "\nCreating reference free alignment for ";
	for (auto it = args.begin(); it != args.end(); ++it){
		cout << *it << " ";
	}
	cout << "\n";
	cout << "Variable sites present in more than " << minproportion*100 << "% of isolates will be included\n";
	cout << "Output will be written to " << outfile << "\n\n";

	alignKmers(minproportion, outfile, args);

	return 0;
}

int compareHelp(void){
	cout << "\nUsage:\n";
	cout << "ska compare [options] <subject split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help\n";
	cout << "-q\t\tQuery split kmer file\n\n";
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

int fastaHelp(void){
	cout << "\nUsage:\n";
	cout << "ska fasta [options] <fasta files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help\n";
	cout << "-k <float>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\t[Default = 15]\n";
	cout << "-o\t\tOutput file name [Default = fasta.kmers]\n\n";
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
					cout << "\nExpecting positive integer after -f flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting positive integer after -f flag\n\n";
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
				outfile = argv[i];
			}
			else {
				cout << "\nExpecting file name after -o flag\n\n";
				return 0;
			}
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

	cout << "\nCreating split " << kmersize << "mers for ";
	for (auto it = args.begin(); it != args.end(); ++it){
		cout << *it << " ";
	}
	cout << "\n";
	cout << "Output will be written to " << outfile << "\n";

	cout << "\n";

	fastaToKmers(args, outfile, kmersize);

	return 0;
}



int fastqHelp(void){
	cout << "\nUsage:\n";
	cout << "ska fastq [options] <fastq files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help\n";
	cout << "-c <float>\tCoverage cutoff. Kmers with coverage below this value will \n\t\tbe discarded. [Default = 8]\n";
	cout << "-f <float>\tFile coverage cutoff. Kmers with coverage below this value \n\t\tin any of the fastq files will be discarded. [Default = 4]\n";
	cout << "-k <float>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\t[Default = 15]\n";
	cout << "-m <float>\tMinimum allowable minor allele frequency. Kmer alleles below \n\t\tthis frequency will be discarded [Default = 0.2]\n";
	cout << "-o\t\tOutput file name [Default = fastq.kmers]\n";
	cout << "-q <float>\tQuality filter for fastq files. No kmers will be created \n\t\tfrom sequence including quality scores below this cutoff. \n\t\t[Default = 30]\n\n";
	return 0;
}


int fastqSubcommand(int argc, char *argv[]){

	string outfile="fastq.kmers";
	long int covcutoff=8;
	long int filecovcutoff=4;
	long int kmersize=15;
	long int minquality=30;
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
		else if (arg == "-f"){
			i++;
			if (i<argc){
				try {
					filecovcutoff = getint(argv[i]);
				}
				catch (const invalid_argument& e) {
					cout << "\nExpecting positive integer after -f flag\n\n";
					return 0;
				}
			}
			else {
				cout << "\nExpecting positive integer after -f flag\n\n";
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
				outfile = argv[i];
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

	cout << "\nCreating split " << kmersize << "mers for ";
	for (auto it = args.begin(); it != args.end(); ++it){
		cout << *it << " ";
	}
	cout << "\n";
	cout << "Filters:\n";
	cout << "\tCoverage cutoff = " << covcutoff << "\n";
	cout << "\tPer file coverage cutoff = " << filecovcutoff << "\n";
	cout << "\tMinimum minor allele frequency = " << minmaf << "\n";
	cout << "\tMinimum base quality = " << minquality << "\n";
	cout << "Output will be written to " << outfile << "\n";

	cout << "\n";

	fastqToKmers(args, outfile, kmersize, minquality, filecovcutoff, covcutoff, minmaf);

	return 0;
}

int mapHelp(void){
	cout << "\nUsage:\n";
	cout << "ska map [options] <split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help\n";
	cout << "-k <float>\tSplit Kmer size. The kmer used for searches will be twice \n\t\tthis length, with the variable base in the middle. e.g. a \n\t\tkmer of 15 will search for 31 base matches with the middle \n\t\tbase being allowed to vary. Must be divisible by 3. \n\t\tMust be the same value used to create the kmer files [Default = 15]\n";
	cout << "-i\t\tInclude reference sequence in alignment\n";
	cout << "-m\t\tMap bases to repeats rather than making them N\n";
	cout << "-o\t\tOutput file name [Default = mappedkmers.aln]\n";
	cout << "-r\t\tReference fasta file name [Required]\n\n";
	return 0;
}


int mapSubcommand(int argc, char *argv[]){

	string outfile = "mappedkmers.aln";
	string reference = "";
	long int kmersize=15;
	bool includeref=false;
	bool maprepeats=false;
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
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = argv[i];
			}
			else {
				cout << "\nExpecting file name after -o flag\n\n";
				return 0;
			}
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

	if (reference==""){
		cout << "\nReference file is required\n";
		mapHelp();
		return 0;
	}

	cout << "\nAligning ";
	for (auto it = args.begin(); it != args.end(); ++it){
		cout << *it << " ";
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

	alignKmersToReference(reference, outfile, args, kmersize, includeref, maprepeats);

	return 0;
}

int summaryHelp(void){
	cout << "\nUsage:\n";
	cout << "ska summary [options] <split kmer files>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help\n\n";
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

int weedHelp(void){
	cout << "\nUsage:\n";
	cout << "ska weed [options] <split kmer file>\n\n";
	cout << "Options:\n";
	cout << "-h\t\tPrint this help\n";
	cout << "-o\t\tOutput file name [Default = weeded.kmers]\n";
	cout << "-w\t\tName of kmer file containing kmers to weed [Required]\n\n";
	return 0;
}


int weedSubcommand(int argc, char *argv[]){

	string outfile="weeded.kmers";
	string weedfile="";
	string kmerfile="";
	vector<string> args;

	for (int i=2; i<argc; ++i){

		string arg=(argv[i]);

		if (arg=="-h" || arg=="--help"){
			weedHelp();
			return 0;
		}
		else if (arg == "-o"){
			i++;
			if (i<argc){
				outfile = argv[i];
			}
			else {
				cout << "\nExpecting file name after -o flag\n\n";
				return 0;
			}
		}
		else if (arg == "-w"){
			i++;
			if (i<argc){
				weedfile = argv[i];
			}
			else {
				cout << "\nExpecting file name after -w flag\n\n";
				return 0;
			}
		}
		else {
			args.push_back(arg);
		}
	}

	if (args.size()!=1){
		cout << "\nExpecting one argument\n";
		weedHelp();
		return 0;
	}

	if (weedfile==""){
		cout << "File containing kmers to weed is required\n";
		weedHelp();
		return 0;
	}

	cout << "\nWeeding " << args[0];
	cout << " to remove kmers in " << weedfile << "\n";
	cout << "Output will be written to " << outfile << "\n";

	cout << "\n";

	weedKmers(args[0], weedfile, outfile);

	return 0;
}


int main(int argc, char *argv[])
{
	
	if (argc<2){
		skaHelp();
		return 0;
	}

	string subcommand=(argv[1]);

	if (subcommand == "align"){
			alignSubcommand(argc, argv);
		}
	else if (subcommand == "compare"){
			compareSubcommand(argc, argv);
		}
	else if (subcommand == "fasta"){
			fastaSubcommand(argc, argv);
		}
	else if (subcommand == "fastq"){
			fastqSubcommand(argc, argv);
		}
	else if (subcommand == "map"){
			mapSubcommand(argc, argv);
		}
	else if (subcommand == "summary"){
			summarySubcommand(argc, argv);
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
	else {
		cout << "\nUnrecognised subcommand\n\n";
		skaHelp();
	}
	
	return 0;

}
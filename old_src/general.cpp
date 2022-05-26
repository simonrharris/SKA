#include <string> //std::string
#include <fstream> //std::ifstream
#include <iostream> //std::cout std::cerr
#include <sstream> //std::stringstream
#include <vector> //std::vector
#include <chrono> //std::chrono
#include "gzstream.h"

std::string versionNumber = "1.0";

std::string splitFileName(const std::string & str){
	std::string filename=str;
	std::string path="";
	std::string::size_type found=filename.find_last_of("/\\");
	if (found!=std::string::npos){
		path=str.substr(0,found);
		filename=str.substr(found+1);
	}
	return filename;
}

int openFileStream(const std::string & fileName, std::ifstream & fileStream, bool verbose){

	fileStream.close();
	fileStream.clear();

	if (verbose){
		std::cout << "Reading " << fileName << std::endl;
	}

	fileStream.open(fileName, std::ios::in);
	if (fileStream.fail()) {
		std::cerr << "Failed to open " << fileName << std::endl << std::endl;
		return 1;
	}
	return 0;
	
}

int openGzFileStream(const std::string & fileName, igzstream & gzFileStream, bool verbose){

	gzFileStream.close();
	gzFileStream.clear();

	gzFileStream.open(fileName.c_str());
	if (gzFileStream.fail()) {
		std::cerr << "Failed to open " << fileName << std::endl << std::endl;
		return 1;
	}
	return 0;
}

int fileToVector(const std::string & filename, std::vector < std::string > & fileargs){
	std::ifstream fileStream;
	fileStream.open(filename, std::ios::in);
	if (fileStream.fail()) {
		std::cerr << "Failed to open " << filename << std::endl << std::endl;
		return 1;
	}
	std::string word;
	while (fileStream >> word){
		if (word.length()>500){
			std::cerr << "Names > 500 characters are not allowed" << std::endl;
			return 1;
		}
		else if (word.length()>0){
			fileargs.push_back(word);
		}
	}
	return 0;
}

void printDuration(const std::chrono::high_resolution_clock::time_point start){
	std::chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Done" << std::endl;
	std::cout << "Total time required: " << elapsed.count() << "s" << std::endl;
}

int calculateMissedSNPs(float totalSiteCount, const float variantSiteCount, const int kmerSize){

	if (variantSiteCount==0){
		std::cout << "Expected number of missed alignments = 0" << std::endl;
		return 0;
	}
	if (totalSiteCount<variantSiteCount){
		std::cout << "Error in calculateMissedSNPs. More variants than sites." << std::endl;
		return 1;
	}

	double nonvariantSiteCount=totalSiteCount-variantSiteCount;
	totalSiteCount--;//as we are assuming the kmer includes one snp
	double snpprobability=nonvariantSiteCount/totalSiteCount;

	for (int i=0; i<kmerSize*2; ++i){
		totalSiteCount--;
		nonvariantSiteCount--;
		snpprobability*=nonvariantSiteCount/totalSiteCount;
	}

	std::cout << "Expected number of missed alignments = " << (1.0-snpprobability)*variantSiteCount << std::endl;
	return 0;
}
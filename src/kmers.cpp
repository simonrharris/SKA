#include <unordered_map> //std::unordered_map
#include <string> //std::string
#include <array> //std::array
#include <vector> // std::vector
#include <sstream> // std::istringstream
#include <cassert> //std::assert
#include <fstream> //std::ifstream
#include <iostream> //std::cout
#include <cmath> //std::pow
#include <set> //std::set
#include "DNA.hpp"
#include "general.hpp"


void ascii_bitstring(std::string & mybits){
	int myremainder=::fmod(int(mybits.length()),6);
	if (myremainder>0){
		for (int i = 0; i<(6-myremainder); ++i){
			mybits.push_back('0');
		}
	}
	for (std::string::size_type i = 0; i<mybits.length(); i+=6){
		assert((i+5)<mybits.length());
		mybits[i/6] = ((int(mybits[i])-'0')+((int(mybits[i+1])-'0')*2)+((int(mybits[i+2])-'0')*4)+((int(mybits[i+3])-'0')*8)+((int(mybits[i+4])-'0')*16)+((int(mybits[i+5])-'0')*32))+33;
	}
	mybits.erase(mybits.length()/6);
}

void vectorbool_from_ascii(std::string & myascii, std::vector < bool > & mybits){
	mybits.resize(myascii.length()*6, false);
	for (std::string::size_type i = 0; i<myascii.length(); ++i){
		int bits=int(myascii[i])-33;
		for (int j=5; j>=0; --j){
			int power=int(pow(2,j));
			mybits[j+(6*i)]=int(bits/power);
			bits=bits%power;
		}
	}
	return;
}

int extractMiddleBase(std::string & myString, char & myChar){
	float stringlen=(myString.length()-1)/2;
	float myremainder=::fmod((myString.length()-1),stringlen);
	if (myremainder>0){
		return 1;
	}
	myChar=myString[stringlen];
	myString.erase(myString.begin()+stringlen);
	return 0;
}

int printKmerFile(const std::unordered_map < std::string, std::array < int, 8 > > & mymap, const std::string outputfile, const int kmersize, float userminmaf){
	
	std::cout << "Writing split kmers to " << outputfile << std::endl;

	std::ofstream kmerfile(outputfile);
	if (kmerfile.fail()){
		std::cerr << std::endl << "Error: Failed to open " << outputfile << std::endl << std::endl;
		return 1;
	}
	kmerfile << "SKA v" << versionNumber << std::endl;
	kmerfile << kmersize << std::endl;
	std:: string samplename=splitFileName(outputfile);
	kmerfile << samplename.substr(0, samplename.find_last_of(".")) << std::endl;
	kmerfile << '"';
	for (std::unordered_map < std::string, std::array < int, 8 > >::const_iterator it=mymap.begin(); it!=mymap.end(); ++it){
		std::string kmer=it->first;
		char base;
		bool basefound=false;
		int i=0;

		int filebasecoverage=0;
		int maxfilebasecoverage=0;
		char maxbase;
		for (std::array < int, 8 >::const_iterator it2=it->second.begin(); it2!=it->second.begin()+4; ++it2, ++i){
			filebasecoverage+=*it2;
			if (*it2>maxfilebasecoverage){
				maxfilebasecoverage=*it2;
				base=bases[i];
			}
		}

		float covcutoff=(1.0-userminmaf)*filebasecoverage;

		if (ascii_codons(kmer)){return 1;}
		if (maxfilebasecoverage>=covcutoff){	
			kmerfile << base << kmer;
		}
		else {
			kmerfile << "N" << kmer;
		}
	}
	kmerfile.close();
	return 0;
}

int printKmerAlleleFrequencies(const std::unordered_map < std::string, std::array < int, 8 > > & mymap, const std::string outputfile){
	
	std::cout << "Writing kmer allele frequencies to " << outputfile << std::endl;

	std::ofstream allelefile(outputfile);
	if (allelefile.fail()){
		std::cerr << std::endl << "Error: Failed to open " << outputfile << std::endl << std::endl;
		return 1;
	}
	for (std::unordered_map < std::string, std::array < int, 8 > >::const_iterator it=mymap.begin(); it!=mymap.end(); ++it){
		std::string kmer=it->first;
		allelefile << kmer;
		for (std::array < int, 8 >::const_iterator it2=it->second.begin(); it2!=it->second.end()-4; ++it2){
		
			allelefile << "\t" << *it2;
		
		}
		allelefile << std::endl;
	}
	allelefile.close();
	return 0;
}

int printMergedKmerFile(const std::unordered_map < std::vector < bool >, std::vector < std::string > > & mymap, const std::string outputfile, const std::vector < std::string > & mysamples, const int kmersize){

	std::cout << "Writing split kmers to " << outputfile << std::endl;

	std::ofstream kmerout(outputfile); //open output file stream
	if (kmerout.fail()){
		std::cerr << std::endl << "Error: Failed to open " << outputfile << std::endl << std::endl;
		return 1;
	}
	kmerout << "SKA v" << versionNumber << std::endl;
	kmerout << kmersize << std::endl; // print kmer size to output file stream
	for ( std::vector < std::string >::const_iterator it=mysamples.begin(); it!=mysamples.end(); ++it){ //print each sample name to output file stream
		kmerout << *it << " "; 
	}
	kmerout << std::endl;

	for ( std::unordered_map < std::vector < bool >, std::vector < std::string > >::const_iterator it=mymap.begin(); it!=mymap.end(); ++it){
		std::stringstream bitstringstream;
		for (std::vector < bool >::const_iterator it2=it->first.begin(); it2!=it->first.end(); ++it2){
			bitstringstream << *it2;
		}
		std::string bitstring = bitstringstream.str();

		ascii_bitstring(bitstring);
		kmerout << bitstring;
		for (std::vector < std::string >::const_iterator it2=it->second.begin(); it2!=it->second.end(); ++it2){
			kmerout << *it2;
		}
		kmerout << std::endl;
	}
	
	kmerout.close();
	return 0;
}


int readKmerHeader(std::ifstream & fileStream, int & kmersize, std::vector < std::string > & names){
	
	std::string line;
	std::getline(fileStream, line);
	fileStream >> kmersize;
	if (kmersize==0){
		return 1;
	}

	std::string name;
	std::getline(fileStream, line);
	std::getline(fileStream, line);

	std::istringstream ss(line);
	while(ss >> name){
		names.push_back(name);
	}

	return 0;
}


int collectSampleNames(const std::vector < std::string > & files, std::vector < std::string > & names, bool verbose){
	if (verbose){
		std::cout << "Collecting number of samples " << std::endl;
	}
	int kmersize;
	std::ifstream fileStream;
	for (std::string::size_type s = 0; s < files.size(); ++s){
		fileStream.open(files[s], std::ios::in);

		if (fileStream.fail()) {
			std::cerr << "Failed to open " << files[s] << std::endl << std::endl;
			return 1;
		}
		int newkmersize;

		if (readKmerHeader(fileStream, newkmersize, names)){return 1;}
		
		if (s==0){
			kmersize=newkmersize;
		}

		if (newkmersize!=kmersize){
			std::cerr << "kmer files have different kmer sizes" << std::endl << std::endl;
			return 1;
		}
		fileStream.close();
	}

	int numSamples=names.size();
	if (verbose){
		std::cout << "Found " << numSamples << " samples in " << files.size() << " files" << std::endl;
	}
	return 0;
}



int getSubsample(const std::vector < std::string > & sample, const std::vector < std::string > & names, std::vector < bool > & include){
	std::set < std::string > sampleSet;

	if (sample.size()>0){
		for (std::vector < std::string >::const_iterator it=sample.begin(); it!=sample.end(); ++it){
			std::pair < std::set < std::string >::iterator,bool> ret = sampleSet.insert(*it);
			if (not ret.second){
				std::cerr << "Warning " << *it << " is in your sample file more than once" << std::endl;
			}
		}
		std::cout << sampleSet.size() << " unique sample names found in your sample file" << std::endl;
	}
	else {
		for (std::vector < std::string >::const_iterator it=names.begin(); it!=names.end(); ++it){
			sampleSet.insert(*it);
		}
	}

	std::vector < bool > included(sampleSet.size());

	for (std::vector < std::string >::const_iterator it=names.begin(); it!=names.end(); ++it){
		std::set < std::string >::iterator it2 = sampleSet.find(*it);
		if (it2 != sampleSet.end()){
			if (included[distance(sampleSet.begin(), it2)]==true){
				std::cerr << "Warning " << *it << " is in your input files more than once. Only the first occurrence will be kept" << std::endl;
			}
			else{
				included[distance(sampleSet.begin(), it2)]=true;
				include[distance(names.begin(), it)]=true;		
			}
		}
	}

	for (std::vector < bool >::iterator it=included.begin(); it!=included.end(); ++it){
		if (not *it){
			std::set < std::string >::iterator it2 = sampleSet.begin();
			advance(it2, distance(included.begin(), it));
			std::cerr << "Warning " << *it2 << " is in your sample file but not in any of your input files" << std::endl;
		}
	}
	std::cout << std::endl;

	return 0;
}


void addKmerToStringMap(std::unordered_map < std::string, std::string > & myKmerMap, const std::string & myKmer, const char myBase, const std::vector < bool > & myBits, const std::vector < int > & mySamples, int currentSampleNumber, int totalSamples, int maxMissing){
	std::unordered_map < std::string, std::string >::iterator kmit = myKmerMap.find(myKmer);//check if the kmer is in the map
	if ( kmit != myKmerMap.end() ){//if the kmer is in the map
		for (int i=0; i<mySamples.size(); ++i){ //add the base to all samples that are true in the bitset
			if (myBits[mySamples[i]]){
				kmit->second[i+currentSampleNumber]=myBase;
			}
		}
	}
	else {//if the kmer isn't in the map
		if ((currentSampleNumber)<=maxMissing){
			std::pair < std::unordered_map < std::string, std::string >::iterator, bool > ret = myKmerMap.insert(std::make_pair(myKmer, std::string (totalSamples,'-')));
			for (int i=0; i<mySamples.size(); ++i){
				if (myBits[mySamples[i]]){
					ret.first->second[i+currentSampleNumber]=myBase;
				}
			}
		}
	}
}


void addKmerToBaseArrayMap(std::unordered_map < std::string, std::array < int, 8 > > & kmerMap, const std::string & kmer, const char base, const bool firstFile){

	std::unordered_map < std::string, std::array < int, 8 > >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
	if ( it != kmerMap.end() ){//if the kmer is in the hash
		it->second[base_score[int(base)]]++;//increment the count of the base for the kmer
	}
	else if (firstFile) {//if the kmers isn't in the hash and we are adding from the first fastq file
		std::pair < std::unordered_map < std::string, std::array < int, 8 > >::iterator,bool>  ret = kmerMap.insert(std::make_pair(kmer, std::array < int, 8 > ()));//insert the kmer into the hash
		ret.first->second[base_score[int(base)]]=1;//increment the count of the base for the kmer
	}
	return;
}


void applyFileKmerArrayMapFilters(std::unordered_map < std::string, std::array < int, 8 > > & kmerMap, const int & userfilecutoff){

	std::cout << "Filtering kmers for file coverage" << std::endl;

	int filebasecoverage;
	int maxfilebasecoverage;
	float filecovcutoff;

	std::unordered_map < std::string, std::array < int, 8 > >::iterator it = kmerMap.begin();
	std::unordered_map < std::string, std::array < int, 8 > >::iterator endIter = kmerMap.end();
	for (; it!=endIter; ){

		//calculate the total file kmer coverage and remove the kmer if it is below the minimum cutoff
		filebasecoverage=0;
		maxfilebasecoverage=0;
		for (std::array < int, 8 >::iterator it2=it->second.begin(); it2!=it->second.begin()+4; ++it2){
			filebasecoverage+=*it2;
			if (*it2>maxfilebasecoverage){
				maxfilebasecoverage=*it2;
			}
		}

		if (filebasecoverage<userfilecutoff){
			kmerMap.erase(it++);
		}
		else {
			int j=4;
			for (int i=0; i<4; ++i, ++j){
				it->second[j]+=it->second[i];
				it->second[i]=0;
			}
			++it;
		}

	}
	std::cout << kmerMap.size() << " unique kmers in map after filtering" << std::endl;
	//return 0;
}


void applyFinalKmerArrayMapFilters(std::unordered_map < std::string, std::array < int, 8 > > & kmerMap, const int & usercovcutoff){

	std::cout << "Filtering kmers for total coverage" << std::endl;

	int basecoverage;
	int maxbasecoverage;
	//float covcutoff;
	float totalcoverage=0;

	std::unordered_map < std::string, std::array < int, 8 > >::iterator it = kmerMap.begin();
	std::unordered_map < std::string, std::array < int, 8 > >::iterator endIter = kmerMap.end();
	for (; it!=endIter; ){

		basecoverage=0;
		maxbasecoverage=0;
		for (std::array < int, 8 >::iterator it2=it->second.begin()+4; it2!=it->second.end(); ++it2){
			basecoverage+=*it2;
			if (*it2>maxbasecoverage){
				maxbasecoverage=*it2;
			}
		}

		//remove kmers with coverage lower than the user defined cutoff
		if (basecoverage<usercovcutoff){
			kmerMap.erase(it++);
		}
		else{
			//filter bases that don't meet the coverage cutoff - don't do this any more
			int j=0;
			for (int i=4; i<8; ++i, ++j){
				/*if (it->second[i]>=covcutoff){
					it->second[j]=it->second[i];
				}*/
				it->second[j]=it->second[i];
				it->second[i]=0;
			}
			totalcoverage+=basecoverage;
			++it;
		}
		
	}

	float meancoverage;
	if (kmerMap.size()>0){
		meancoverage=totalcoverage/kmerMap.size();
	}
	else{
		meancoverage=0;
	}
	std::cout << kmerMap.size() << " unique kmers in map after filtering" << std::endl;
	std::cout << "Mean kmer coverage is " << meancoverage << std::endl;

	//return 0;
}


void reverseVectorBoolKmerMap(std::unordered_map < std::string, std::vector < bool > > & kmerMap, std::unordered_map < std::vector < bool >, std::vector < std::string > > & revKmerMap){

	std::unordered_map < std::string, std::vector < bool > >::iterator it = kmerMap.begin();
	std::unordered_map < std::string, std::vector < bool > >::iterator endIter = kmerMap.end();

	for (; it!=endIter; ){
		std::unordered_map < std::vector < bool >, std::vector < std::string > >::iterator it2 = revKmerMap.find(it->second);//check if the bitset is in the map
		if ( it2 != revKmerMap.end() ){//if the bitset is in the map
			it2->second.push_back(it->first); //add the kmer
		}
		else {//if the bitset isn't in the map
			revKmerMap.insert(std::make_pair(it->second, std::vector < std::string > {it->first})); //insert the vector into the map
		}
		kmerMap.erase(it++); //remove kmer from kmerMap to save space
	}

}


void reverseStringKmerMap(std::unordered_map < std::string, std::string > & kmerMap, std::unordered_map < std::vector < bool >, std::vector < std::string > > & revKmerMap){

	std::unordered_map < std::string, std::string >::iterator it = kmerMap.begin();
	std::unordered_map < std::string, std::string >::iterator endIter = kmerMap.end();

	for (; it!=endIter; ){

		for (int i=0; i<5; ++i){
			std::string kmer (1, bases[i]);
			kmer+=it->first;
			std::vector < bool > myvb (it->second.length(), false);

			int j=0;
			bool matched=false;
			for (std::string::iterator sit=it->second.begin(); sit!=it->second.end(); ++sit, ++j){
				if (*sit==bases[i]){
					myvb[j]=true;
					matched=true;
				}
			}

			if (matched){
				std::unordered_map < std::vector < bool >, std::vector < std::string > >::iterator it2 = revKmerMap.find(myvb);//check if the bitset is in the map
				if ( it2 != revKmerMap.end() ){//if the bitset is in the map
					it2->second.push_back(kmer); //add the kmer
				}
				else {//if the bitset isn't in the map
					revKmerMap.insert(std::make_pair(myvb, std::vector < std::string > {kmer})); //insert the vector into the map
				}
			}
		}
		kmerMap.erase(it++); //remove kmer from kmerMap to save space
	}

}

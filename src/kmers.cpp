#include <unordered_map> //std::unordered_map
#include <string> //std::string
#include <array> //std::array
#include <vector> // std::vector
#include <sstream> // std::istringstream
#include <chrono> //timing
#include <cassert> //std::assert
#include <fstream> //std::ifstream
#include <iostream> //std::cout
#include <cmath> //std::pow
#include <set> //std::set
#include "DNA.hpp"

using namespace std;


void ascii_bitstring(string & mybits){
	int myremainder=::fmod(int(mybits.length()),6);
	if (myremainder>0){
		for (int i = 0; i<(6-myremainder); ++i){
			mybits.push_back('0');
		}
	}
	for (string::size_type i = 0; i<mybits.length(); i+=6){
		assert((i+5)<mybits.length());
		mybits[i/6] = ((int(mybits[i])-'0')+((int(mybits[i+1])-'0')*2)+((int(mybits[i+2])-'0')*4)+((int(mybits[i+3])-'0')*8)+((int(mybits[i+4])-'0')*16)+((int(mybits[i+5])-'0')*32))+33;
	}
	mybits.erase(mybits.length()/6);
}

void vectorbool_from_ascii(string & myascii, vector < bool > & mybits){
	mybits.resize(myascii.length()*6, false);
	for (string::size_type i = 0; i<myascii.length(); ++i){
		int bits=int(myascii[i])-33;
		for (int j=5; j>=0; --j){
			int power=int(pow(2,j));
			mybits[j+(6*i)]=int(bits/power);
			bits=bits%power;
		}
	}
	return;
}

int extractMiddleBase(string & myString, char & myChar){
	float stringlen=(myString.length()-1)/2;
	float myremainder=::fmod((myString.length()-1),stringlen);
	if (myremainder>0){
		return 1;
	}
	myChar=myString[stringlen];
	myString.erase(myString.begin()+stringlen);
	return 0;
}

int printKmerFile(const unordered_map <string, array < int, 8 > > & mymap, const string outputfile, const int kmersize){
	
	cout << "Writing kmers to " << outputfile << endl;

	ofstream kmerfile(outputfile);
	if (kmerfile.fail()){
		cerr << endl << "Error: Failed to open " << outputfile << endl << endl;
		return 1;
	}
	
	kmerfile << kmersize << endl;
	kmerfile << outputfile.substr(0, outputfile.find_last_of(".")) << endl;
	kmerfile << '"';
	for (unordered_map <string, array < int, 8 > >::const_iterator it=mymap.begin(); it!=mymap.end(); ++it){
		string kmer=it->first;
		char base;
		bool basefound=false;
		int i=0;
		for (array < int, 8 >::const_iterator it2=it->second.begin(); it2!=it->second.end()-4; ++it2, ++i){
		
			if (*it2>0){
				if (basefound){
					base='N';
					break;
				}
				else{
					base=bases[i];
					basefound=true;
				}
			}
		
		}
		if (basefound){
			ascii_codons(kmer);
			kmerfile << base << kmer;
		}
		else {
			cout << "Error: Failed to print kmer with no middle base" << endl;
			return 1;
		}
	}
	kmerfile.close();
	return 0;
}

int printMergedKmerFile(const unordered_map < vector < bool >, vector < string > > & mymap, const string outfileprefix, const vector < string > & mysamples, const int kmersize){

	string outputfile;

	if (mysamples.size()==1){
		outputfile=outfileprefix+".kmers";
	}
	else {
		outputfile=outfileprefix+".kmerge";
	}

	cout << "Writing merged file to " << outputfile << endl;

	ofstream kmerout(outputfile); //open output file stream
	if (kmerout.fail()){
		cerr << endl << "Error: Failed to open " << outputfile << endl << endl;
		return 1;
	}

	kmerout << kmersize << endl; // print kmer size to output file stream
	for ( vector < string >::const_iterator it=mysamples.begin(); it!=mysamples.end(); ++it){ //print each sample name to output file stream
		kmerout << *it << " "; 
	}
	kmerout << endl;

	for ( unordered_map < vector < bool >, vector < string > >::const_iterator it=mymap.begin(); it!=mymap.end(); ++it){
		stringstream bitstringstream;
		for (vector < bool >::const_iterator it2=it->first.begin(); it2!=it->first.end(); ++it2){
			bitstringstream << *it2;
		}
		string bitstring = bitstringstream.str();

		ascii_bitstring(bitstring);
		kmerout << bitstring;
		for (vector < string >::const_iterator it2=it->second.begin(); it2!=it->second.end(); ++it2){
			kmerout << *it2;
		}
		kmerout << endl;
	}
	
	kmerout.close();
	return 0;
}


int readKmerHeader(ifstream & fileStream, int & kmersize, vector < string > & names){
	
	string line;
	fileStream >> kmersize;
	if (kmersize==0){
		return 1;
	}

	string name;
	getline(fileStream, line);
	getline(fileStream, line);

	istringstream ss(line);
	while(ss >> name){
		names.push_back(name);
	}

	return 0;
}


int collectSampleNames(const vector < string > & files, vector < string > & names){
	cout << "Collecting number of samples " << endl;
	int kmersize;
	ifstream fileStream;
	for (string::size_type s = 0; s < files.size(); ++s){
		fileStream.open(files[s], ios::in);

		if (fileStream.fail()) {
			cerr << "Failed to open " << files[s] << endl << endl;
			return 1;
		}
		int newkmersize;

		if (readKmerHeader(fileStream, newkmersize, names)){return 1;}
		
		if (s==0){
			kmersize=newkmersize;
		}

		if (newkmersize!=kmersize){
			cerr << "kmer files have different kmer sizes" << endl << endl;
			return 1;
		}
		fileStream.close();
	}

	int numSamples=names.size();

	cout << "Found " << numSamples << " samples in " << files.size() << " files" << endl;
	return 0;
}



int getSubsample(const vector < string > & sample, const vector < string > & names, vector < bool > & include){
	set < string > sampleSet;

	if (sample.size()>0){
		for (vector < string >::const_iterator it=sample.begin(); it!=sample.end(); ++it){
			pair<set<string>::iterator,bool> ret = sampleSet.insert(*it);
			if (not ret.second){
				cerr << "Warning " << *it << " is in your sample file more than once" << endl;
			}
		}
		cout << sampleSet.size() << " unique sample names found in your sample file" << endl;
	}
	else {
		for (vector < string >::const_iterator it=names.begin(); it!=names.end(); ++it){
			sampleSet.insert(*it);
		}
	}

	vector < bool > included(sampleSet.size());

	for (vector < string >::const_iterator it=names.begin(); it!=names.end(); ++it){
		set < string >::iterator it2 = sampleSet.find(*it);
		if (it2 != sampleSet.end()){
			if (included[distance(sampleSet.begin(), it2)]==true){
				cerr << "Warning " << *it << " is in your input files more than once. Only the first occurrence will be kept" << endl;
			}
			else{
				included[distance(sampleSet.begin(), it2)]=true;
				include[distance(names.begin(), it)]=true;		
			}
		}
	}

	for (vector < bool >::iterator it=included.begin(); it!=included.end(); ++it){
		if (not *it){
			set < string >::iterator it2 = sampleSet.begin();
			advance(it2, distance(included.begin(), it));
			cerr << "Warning " << *it2 << " is in your sample file but not in any of your input files" << endl;
		}
	}
	cout << endl;

	return 0;
}


void addKmerToBaseArrayMap(unordered_map <string, array < int, 8 > > & kmerMap, const string & kmer, const char base, const bool firstFile){

	unordered_map <string, array < int, 8 > >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
	if ( it != kmerMap.end() ){//if the kmer is in the hash
		it->second[base_score[base]]++;//increment the count of the base for the kmer
	}
	else if (firstFile) {//if the kmers isn't in the hash and we are adding from the first fastq file
		pair<unordered_map <string, array < int, 8 > >::iterator,bool>  ret = kmerMap.insert(make_pair(kmer, array < int, 8 > ()));//insert the kmer into the hash
		ret.first->second[base_score[base]]=1;//increment the count of the base for the kmer
	}
	return;
}


int applyFileKmerArrayMapFilters(unordered_map <string, array < int, 8 > > & kmerMap, const int & userfilecutoff, const float & userminmaf){

	cout << "Filtering kmers for file coverage" << endl;

	int filebasecoverage;
	int maxfilebasecoverage;
	float filecovcutoff;

	unordered_map <string, array < int, 8 > >::iterator it = kmerMap.begin();
	unordered_map <string, array < int, 8 > >::iterator endIter = kmerMap.end();
	for (; it!=endIter; ){

		//calculate the total file kmer coverage and remove the kmer if it is below the minimum cutoff
		filebasecoverage=0;
		maxfilebasecoverage=0;
		for (array < int, 8 >::iterator it2=it->second.begin(); it2!=it->second.begin()+4; ++it2){
			filebasecoverage+=*it2;
			if (*it2>maxfilebasecoverage){
				maxfilebasecoverage=*it2;
			}
		}

		filecovcutoff=float(filebasecoverage)*userminmaf;

		if (maxfilebasecoverage<filecovcutoff || filebasecoverage<userfilecutoff){
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
	cout << kmerMap.size() << " unique kmers in map after filtering" << endl;
	return 0;
}


int applyFinalKmerArrayMapFilters(unordered_map <string, array < int, 8 > > & kmerMap, const int & usercovcutoff, const float & userminmaf){

	cout << "Filtering kmers for total coverage" << endl;

	int basecoverage;
	int maxbasecoverage;
	float covcutoff;
	float totalcoverage=0;

	unordered_map <string, array < int, 8 > >::iterator it = kmerMap.begin();
	unordered_map <string, array < int, 8 > >::iterator endIter = kmerMap.end();
	for (; it!=endIter; ){

		basecoverage=0;
		maxbasecoverage=0;
		for (array < int, 8 >::iterator it2=it->second.begin()+4; it2!=it->second.end(); ++it2){
			basecoverage+=*it2;
			if (*it2>maxbasecoverage){
				maxbasecoverage=*it2;
			}
		}

		covcutoff=float(basecoverage)*userminmaf;
		/*if (covcutoff<usercovcutoff){
			covcutoff=usercovcutoff;
		}*/

		//remove kmers with coverage lower than the user defined cutoff
		if (maxbasecoverage<covcutoff || basecoverage<usercovcutoff){
			kmerMap.erase(it++);
		}
		else{
			//filter bases that don't meet the coverage cutoff
			int j=0;
			for (int i=4; i<8; ++i, ++j){
				if (it->second[i]>=covcutoff){
					it->second[j]=it->second[i];
				}
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
	cout << kmerMap.size() << " unique kmers in map after filtering" << endl;
	cout << "Mean kmer coverage is " << meancoverage << endl;

	return 0;
}



void reverseVectorBoolKmerMap(unordered_map < string, vector < bool > > & kmerMap, unordered_map < vector < bool >, vector < string > > & revKmerMap){

	unordered_map < string, vector < bool > >::iterator it = kmerMap.begin();
	unordered_map < string, vector < bool > >::iterator endIter = kmerMap.end();

	for (; it!=endIter; ){
		unordered_map < vector < bool >, vector < string > >::iterator it2 = revKmerMap.find(it->second);//check if the bitset is in the map
		if ( it2 != revKmerMap.end() ){//if the bitset is in the map
			it2->second.push_back(it->first); //add the kmer
		}
		else {//if the bitset isn't in the map
			vector <string> myvector; //create a new vector of strings
			myvector.push_back(it->first); //add the bitset to the vector
			revKmerMap.insert(make_pair(it->second, myvector)); //insert the vector into the map
		}
		kmerMap.erase(it++); //remove kmer from kmerMap to save space
	}

}

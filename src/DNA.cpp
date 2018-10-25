#include <string> //std::string
#include <algorithm> //std::reverse std::transform
#include <sstream> // std::stringstream
#include <iostream> //std:cerr
#include "gzstream.h"


char complement_table[128] = {
      '-',
	  '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
	  '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
	  '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
	  '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
	  'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O','P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', 
	  '-', '-', '-', '-', '-', '-', 
	  't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o','p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 
	  '-', '-', '-', '-', '-'
};

char base_score[128] = {
      5,
	  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	  0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 
	  5, 5, 5, 5, 5, 5, 
	  0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4,
	  5, 5, 5, 5, 5
};

char complement_base_score[5] = {
      4, 3, 2, 1, 5
};

char bases[5] = {
          'A', 'C', 'G', 'T', 'N'
};

char complement_bases[5] = {
          'T', 'G', 'C', 'A', 'N'
};

char complement(const char & n){   
	return complement_table[int(n)];
}

bool reverseComplementIsMin(const std::string & mystring){
	for (float i = 0; i < ((mystring.length())/2); ++i){
		int j = i+1;
		if (mystring[i]<complement_table[int(mystring[mystring.length()-j])]){
			return false;
		}
		else if (complement_table[int(mystring[mystring.length()-j])]<mystring[i]){
			return true;
		}
	}
	return false;
}

bool reverseComplementIfMin(std::string & mystring){
	if (reverseComplementIsMin(mystring)){
		std::reverse(mystring.begin(), mystring.end());
		std::transform(mystring.begin(),mystring.end(),mystring.begin(),complement);
		return true;
	}
	return false;
}

int lowqualitytoN(std::string & mysequence,const std::string & myquality, const int & userminquality, int adjustment){
	int minquality=userminquality+adjustment;
	for (std::string::size_type i = 0; i<mysequence.length(); ++i){
		if (int(myquality[i])<minquality){
			mysequence[i]='N';
		}
	}
	return 0;
}

int IUPACToN(std::string & mysequence){
	for (std::string::iterator it=mysequence.begin(); it!=mysequence.end(); ++it){
		if (base_score[int(*it)]>4){
			std::cerr << "Error: Unrecognised character " << *it << "in sequence" << std::endl << std::endl;
			return 1;
		}
		else if (base_score[int(*it)]>3){
			*it='N';
		}
	}
	return 0;
}

int ascii_codons(std::string & myDNA){
	for (std::string::size_type i = 0; i<myDNA.length(); i+=3){
		if ((i+2)>=myDNA.length()){
			std::cerr << "Kmer length must be divisible by 3" << std::endl << std::endl;
			return 1;
		}
		myDNA[i/3] = (base_score[int(myDNA[i])]+(base_score[int(myDNA[i+1])]*4)+(base_score[int(myDNA[i+2])]*16))+63;
	}
	myDNA.erase(myDNA.length()/3);
	return 0;
}

std::string codons_from_ascii(std::string & myascii){
	std::string myDNA ((myascii.length()*3), '-');
	for (std::string::size_type i = 0; i<myascii.length(); ++i){
		myDNA[(i*3)+2]=bases[int(myascii[i]-63)/16];
		int sixteensreamainder=int(myascii[i]-63)%16;
		myDNA[(i*3)+1]=bases[sixteensreamainder/4];
		myDNA[(i*3)]=bases[sixteensreamainder%4];
	}
	return myDNA;
}


int readNextFastaSequence(igzstream & gzfileStream, const std::string & filename, std::string & name, std::string & sequence){
	std::string line;
	line.reserve(10000);
	if (gzfileStream.get()!='>'){//check the fist character is a > and remove it
		std::cerr << "Error: " << filename << " is not in the correct format. Expecting header line to start with >." << std::endl << std::endl;
		return 1;
	}
  	std::getline(gzfileStream, line);//read the rest of the headerline
  	std::stringstream headerstream;
	headerstream << line;//convert the sequence to stringstream
	std::getline(headerstream, name, ' ');
  	while (gzfileStream.peek()!=EOF && gzfileStream.peek()!='>'){//read the sequence
  		std::getline(gzfileStream, line);
  		sequence.append(line);
  	}
	if (gzfileStream.fail()){
		std::cerr << filename << " is not in the correct format." << std::endl << std::endl;
		return 1;
	}
	sequence.erase(std::remove_if( sequence.begin(), sequence.end(), ::isspace ), sequence.end() );//remove any gaps from the sequence
	std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);//change all letters in the string to upper case
	return 0;
}


int readNextFastqSequence(igzstream & gzfileStream, const std::string & filename, std::string & sequence, std::string & quality){
	std::string line;
	line.reserve(1000);
	if (gzfileStream.peek()!='@'){
		std::cerr << "Error: " << filename << " is not in the correct format. Expecting header line to start with @." << std::endl << std::endl;
		return 1;
	}
  	getline(gzfileStream, line);
  	getline(gzfileStream, sequence);
  	if (gzfileStream.peek()!='+'){
		std::cerr << "Error: " << filename << " is not in the correct format. Expecting separator line to start with +." << std::endl << std::endl;
		return 1;
	}
  	getline(gzfileStream, line);
  	getline(gzfileStream, quality);
  	if (quality.length()!=sequence.length()){
  		std::cerr << "Error: " << filename << " is not in the correct format. Sequence and quality lines must be of equal length." << std::endl << std::endl;
		return 1;
  	}
	if (gzfileStream.fail()){
		std::cerr << "Error: " << filename << " is not in the correct format. Expecting a multiple of 4 lines." << std::endl << std::endl;
		return 1;
	}
	std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);//change all letters in the string to upper case	
	return 0;
}


int countSequencesinFasta(const std::string & filename, int & sequenceCount){	
	igzstream gzfileStream;
	gzfileStream.open(filename.c_str());
	if (gzfileStream.get()!='>'){
		std::cerr << "Error: " << filename << " is not in the correct format. Expecting header to start with >." << std::endl << std::endl;
		return 1;
	}
	std::string sequence;
	while (getline(gzfileStream, sequence, '>')){
		sequenceCount++;
	}
	gzfileStream.close();
	return 0;
}

int circulariseSequence(std::string & mySequence, const int myKmerLength){

	if (mySequence.length()<myKmerLength*2){
		return 1;
	}

	std::string myNewEnd=mySequence.substr(0,myKmerLength);
	std::string myNewStart=mySequence.substr(mySequence.length()-myKmerLength,myKmerLength);
	mySequence.reserve(mySequence.length()+(myKmerLength*2));
	mySequence.insert(mySequence.length(), myNewEnd);
	mySequence.insert(0, myNewStart);

	return 0;

}


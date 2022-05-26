// Copyright 2005, Google Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "../src/DNA.hpp"
#include "gtest/gtest.h"
#include <string>
namespace {

  // Tests complement

  TEST(complementTest, bases){

    EXPECT_EQ('A', complement('T'));
    EXPECT_EQ('C', complement('G'));
    EXPECT_EQ('G', complement('C'));
    EXPECT_EQ('T', complement('A'));
    EXPECT_EQ('A', complement('U'));
    EXPECT_EQ('a', complement('t'));
    EXPECT_EQ('c', complement('g'));
    EXPECT_EQ('g', complement('c'));
    EXPECT_EQ('t', complement('a'));
    EXPECT_EQ('a', complement('u'));

  }

  TEST(complementTest, IUPAC){

    EXPECT_EQ('R', complement('Y'));
    EXPECT_EQ('Y', complement('R'));
    EXPECT_EQ('S', complement('S'));
    EXPECT_EQ('W', complement('W'));
    EXPECT_EQ('K', complement('M'));
    EXPECT_EQ('M', complement('K'));
    EXPECT_EQ('K', complement('M'));
    EXPECT_EQ('B', complement('V'));
    EXPECT_EQ('D', complement('H'));
    EXPECT_EQ('H', complement('D'));
    EXPECT_EQ('V', complement('B'));
    EXPECT_EQ('r', complement('y'));
    EXPECT_EQ('y', complement('r'));
    EXPECT_EQ('s', complement('s'));
    EXPECT_EQ('w', complement('w'));
    EXPECT_EQ('k', complement('m'));
    EXPECT_EQ('m', complement('k'));
    EXPECT_EQ('b', complement('v'));
    EXPECT_EQ('d', complement('h'));
    EXPECT_EQ('h', complement('d'));
    EXPECT_EQ('v', complement('b'));
    EXPECT_EQ('n', complement('n'));

  }

    // Tests reverseComplementIsMin

  TEST(reverseComplementIsMinTest, reverseIsMin){

    EXPECT_EQ(true, reverseComplementIsMin("GAG"));
    EXPECT_EQ(true, reverseComplementIsMin("GGAGG"));
    EXPECT_EQ(true, reverseComplementIsMin("GGGAGGG"));
    EXPECT_EQ(true, reverseComplementIsMin("GGGGAGGGG"));
    EXPECT_EQ(true, reverseComplementIsMin("GGGGGAGGGGG"));
    EXPECT_EQ(true, reverseComplementIsMin("TAT"));
    EXPECT_EQ(true, reverseComplementIsMin("TTATT"));
    EXPECT_EQ(true, reverseComplementIsMin("TTTATTT"));
    EXPECT_EQ(true, reverseComplementIsMin("TTTTATTTT"));
    EXPECT_EQ(true, reverseComplementIsMin("TTTTTATTTTT"));
    EXPECT_EQ(true, reverseComplementIsMin("AGAGT"));
    EXPECT_EQ(true, reverseComplementIsMin("AAGAGTT"));
    EXPECT_EQ(true, reverseComplementIsMin("AAAGAGTTT"));
    EXPECT_EQ(true, reverseComplementIsMin("AAAAGAGTTTT"));
    EXPECT_EQ(true, reverseComplementIsMin("AAAAAGAGTTTTT"));
    EXPECT_EQ(true, reverseComplementIsMin("CGAGG"));
    EXPECT_EQ(true, reverseComplementIsMin("CCGAGGG"));
    EXPECT_EQ(true, reverseComplementIsMin("CCCGAGGGG"));
    EXPECT_EQ(true, reverseComplementIsMin("CCCCGAGGGGG"));
    EXPECT_EQ(true, reverseComplementIsMin("CCCCCGAGGGGGG"));

  }

  TEST(reverseComplementIsMinTest, reverseIsNotMin){

    EXPECT_EQ(false, reverseComplementIsMin("CAC"));
    EXPECT_EQ(false, reverseComplementIsMin("CCACC"));
    EXPECT_EQ(false, reverseComplementIsMin("CCCACCC"));
    EXPECT_EQ(false, reverseComplementIsMin("CCCCACCCC"));
    EXPECT_EQ(false, reverseComplementIsMin("CCCCCACCCCC"));
    EXPECT_EQ(false, reverseComplementIsMin("AAA"));
    EXPECT_EQ(false, reverseComplementIsMin("AAAAA"));
    EXPECT_EQ(false, reverseComplementIsMin("AAAAAAA"));
    EXPECT_EQ(false, reverseComplementIsMin("AAAAAAAAA"));
    EXPECT_EQ(false, reverseComplementIsMin("AAAAAAAAAAA"));
    EXPECT_EQ(false, reverseComplementIsMin("TCACA"));
    EXPECT_EQ(false, reverseComplementIsMin("TTCACAA"));
    EXPECT_EQ(false, reverseComplementIsMin("TTTCACAAA"));
    EXPECT_EQ(false, reverseComplementIsMin("TTTTCACAAAA"));
    EXPECT_EQ(false, reverseComplementIsMin("TTTTTCACAAAAA"));
    EXPECT_EQ(false, reverseComplementIsMin("GCACC"));
    EXPECT_EQ(false, reverseComplementIsMin("GGCACCC"));
    EXPECT_EQ(false, reverseComplementIsMin("GGGCACCCC"));
    EXPECT_EQ(false, reverseComplementIsMin("GGGGCACCCCC"));
    EXPECT_EQ(false, reverseComplementIsMin("GGGGGCACCCCCC"));

  }

  TEST(reverseComplementIsMinTest, palindrome){

    EXPECT_EQ(false, reverseComplementIsMin("AAATT"));
    EXPECT_EQ(false, reverseComplementIsMin("AACTT"));
    EXPECT_EQ(false, reverseComplementIsMin("AAGTT"));
    EXPECT_EQ(false, reverseComplementIsMin("AATTT"));

  }

    // Tests lowqualitytoN

  TEST(lowqualitytoNTest, belowMinQuality){

    std::string sequence="AAAAAAAAAAAAAAAAAA";
    std::string quality="HHHHHHH4H>HHHHHHHH";
    int minqual=20;
    EXPECT_EQ(0, lowqualitytoN(sequence, quality, minqual));
    EXPECT_STREQ("AAAAAAANAAAAAAAAAA", sequence.c_str());
    minqual=30;
    EXPECT_EQ(0, lowqualitytoN(sequence, quality, minqual));
    EXPECT_STREQ("AAAAAAANANAAAAAAAA", sequence.c_str());
    minqual=40;
    EXPECT_EQ(0, lowqualitytoN(sequence, quality, minqual));
    EXPECT_STREQ("NNNNNNNNNNNNNNNNNN", sequence.c_str());

  }

  TEST(lowqualitytoNTest, aboveMinQuality){

    std::string sequence="AAAAAAAAAAAAAAAAAA";
    std::string quality="IIIIIII5I?IIIIIIII";
    int minqual=20;
    EXPECT_EQ(0, lowqualitytoN(sequence, quality, minqual));
    EXPECT_STREQ("AAAAAAAAAAAAAAAAAA", sequence.c_str());
    minqual=30;
    EXPECT_EQ(0, lowqualitytoN(sequence, quality, minqual));
    EXPECT_STREQ("AAAAAAANAAAAAAAAAA", sequence.c_str());
    minqual=40;
    EXPECT_EQ(0, lowqualitytoN(sequence, quality, minqual));
    EXPECT_STREQ("AAAAAAANANAAAAAAAA", sequence.c_str());

  }

  TEST(lowqualitytoNTest, noLowQuality){

    std::string sequence="AAAAAAAAAAAAAAAAAA";
    std::string quality="ZZZZZZZZZZZZZZZZZZ";
    int minqual=20;
    EXPECT_EQ(0, lowqualitytoN(sequence, quality, minqual));
    EXPECT_STREQ("AAAAAAAAAAAAAAAAAA", sequence.c_str());
    minqual=30;
    EXPECT_EQ(0, lowqualitytoN(sequence, quality, minqual));
    EXPECT_STREQ("AAAAAAAAAAAAAAAAAA", sequence.c_str());
    minqual=40;
    EXPECT_EQ(0, lowqualitytoN(sequence, quality, minqual));
    EXPECT_STREQ("AAAAAAAAAAAAAAAAAA", sequence.c_str());

  }

  // Test IUPACToN

  TEST(IUPACToNTest, allACGT){

    std::string sequence="ACGT";
    EXPECT_EQ(0, IUPACToN(sequence));
    EXPECT_STREQ("ACGT", sequence.c_str());

  }

  TEST(IUPACToNTest, withN){

    std::string sequence="ACNGT";
    EXPECT_EQ(0, IUPACToN(sequence));
    EXPECT_STREQ("ACNGT", sequence.c_str());

  }

  TEST(IUPACToNTest, withIUPAC){

    std::string sequence="ACRYSWKMBDHVNGT";
    EXPECT_EQ(0, IUPACToN(sequence));
    EXPECT_STREQ("ACNNNNNNNNNNNGT", sequence.c_str());

  }

  TEST(IUPACToNTest, nonIUPAC){

    std::string sequence="AC-GT";
    EXPECT_EQ(1, IUPACToN(sequence));
    sequence="AC.GT";
    EXPECT_EQ(1, IUPACToN(sequence));
    sequence="AC;GT";
    EXPECT_EQ(1, IUPACToN(sequence));
    sequence="AC,GT";
    EXPECT_EQ(1, IUPACToN(sequence));
    sequence="AC]GT";
    EXPECT_EQ(1, IUPACToN(sequence));
    sequence="AC^GT";
    EXPECT_EQ(1, IUPACToN(sequence));

  }

  // Test ascii_codons

  TEST(ascii_codonsTest, codonsToAscii){

    std::string sequence="ACTTGCTACTTTAAGCGA";
    ascii_codons(sequence);
    EXPECT_STREQ("sZR~_H", sequence.c_str());

  }



    // Test codons_from_ascii

  TEST(codons_from_asciiTest, codonsFromAscii){

    std::string asciisequence="sZR~_H";
    codons_from_ascii(asciisequence);
    EXPECT_STREQ("sZR~_H", asciisequence.c_str());

  }

  // Test circulariseSequence

  TEST(circulariseSequenceTest, circulariseSequence){

    std::string sequence="AAACCCCCCCCTTT";
    circulariseSequence(sequence, 3);
    EXPECT_STREQ("TTTAAACCCCCCCCTTTAAA", sequence.c_str());

  }

} 

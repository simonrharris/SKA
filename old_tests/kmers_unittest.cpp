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


#include "../src/kmers.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
namespace {

  // Tests ascii_bitstring

  TEST(ascii_bitstringTest, divisibleBy6){

    std::string bitstring= "000000";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("!", bitstring.c_str());
    bitstring= "111111";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("`", bitstring.c_str());
    bitstring= "000000000000";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("!!", bitstring.c_str());
    bitstring= "111111111111";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("``", bitstring.c_str());

  }

  TEST(ascii_bitstringTest, notDivisibleBy6){

    std::string bitstring= "00000000000";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("!!", bitstring.c_str());
    bitstring= "0000000000";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("!!", bitstring.c_str());
    bitstring= "000000000";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("!!", bitstring.c_str());
    bitstring= "00000000";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("!!", bitstring.c_str());
    bitstring= "0000000";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("!!", bitstring.c_str());
    bitstring= "00000";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("!", bitstring.c_str());
    bitstring= "0000";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("!", bitstring.c_str());
    bitstring= "000";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("!", bitstring.c_str());
    bitstring= "00";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("!", bitstring.c_str());
    bitstring= "0";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("!", bitstring.c_str());
    bitstring= "11111";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("@", bitstring.c_str());
    bitstring= "1111";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("0", bitstring.c_str());
    bitstring= "111";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("(", bitstring.c_str());
    bitstring= "11";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("$", bitstring.c_str());
    bitstring= "1";
    ascii_bitstring(bitstring);
    EXPECT_STREQ("\"", bitstring.c_str());

  }


  // Tests vectorbool_from_ascii

  TEST(vectorbool_from_asciiTest, allFalse){

    std::vector < bool > myvectorbool;
    std::string asciistring= "!";
    vectorbool_from_ascii(asciistring, myvectorbool);
    EXPECT_EQ(6, myvectorbool.size());
    for (std::vector < bool >::iterator it=myvectorbool.begin(); it!=myvectorbool.end(); ++it){//There is probably a better way to test euqlity of this
      EXPECT_EQ(false, *it);
    }
    myvectorbool.clear();
    asciistring= "!!";
    vectorbool_from_ascii(asciistring, myvectorbool);
    EXPECT_EQ(12, myvectorbool.size());
    for (std::vector < bool >::iterator it=myvectorbool.begin(); it!=myvectorbool.end(); ++it){
      EXPECT_EQ(false, *it);
    }

  }

    TEST(vectorbool_from_asciiTest, allTrue){

    std::vector < bool > myvectorbool;
    std::string asciistring= "`";
    vectorbool_from_ascii(asciistring, myvectorbool);
    EXPECT_EQ(6, myvectorbool.size());
    for (std::vector < bool >::iterator it=myvectorbool.begin(); it!=myvectorbool.end(); ++it){
      EXPECT_EQ(true, *it);
    }
    myvectorbool.clear();
    asciistring= "``";
    vectorbool_from_ascii(asciistring, myvectorbool);
    EXPECT_EQ(12, myvectorbool.size());
    for (std::vector < bool >::iterator it=myvectorbool.begin(); it!=myvectorbool.end(); ++it){
      EXPECT_EQ(true, *it);
    }

  }

  // Tests extractMiddleBase

  TEST(extractMiddleBaseTest, oddBases){
    char testbase;
    std::string teststring="TCTATCT";
    EXPECT_EQ(0, extractMiddleBase(teststring, testbase));
    EXPECT_EQ('A', testbase);
    teststring="TCTGTCTCT";
    EXPECT_EQ(0, extractMiddleBase(teststring, testbase));
    EXPECT_EQ('T', testbase);

  }

 TEST(extractMiddleBaseTest, evenBases){
    char testbase;
    std::string teststring="TCATCT";
    EXPECT_EQ(1, extractMiddleBase(teststring, testbase));
    teststring="TCAGGTCT";
    EXPECT_EQ(1, extractMiddleBase(teststring, testbase));

  }


} 

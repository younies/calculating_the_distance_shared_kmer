//
//  enumKTree.hpp
//  calculating_distance
//
//  Created by Younies Mahmoud on 5/2/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef enumKTree_hpp
#define enumKTree_hpp

#include <stdio.h>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <memory.h>
#include <fstream>
#include "helpers.hpp"

using namespace std;




struct enumNode {
    pair<int, int>  I;
    pair<int, int> I_dag;
    int depth;
    int base;
} ;

struct BWT2 {
    string BWT;
    string BWT_dag;
    vector<vector<int> > occ_B;
    vector<vector<int> > occ_reverse_B;
    vector<int> C;
};

struct Distance {
    long long count;
    long long first_DNA_count;
    long long second_DNA_count;
} ;


Distance enum_difference(int k , BWT2 &DNA1 , BWT2 &DNA2 );


bool computing_child(enumNode &parent, char base , BWT2 & DNA , enumNode &newChild );

bool has_interval(enumNode node);

void buid_bwt(BWT2 & DNA);


#endif /* enumKTree_hpp */

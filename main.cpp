//
//  main.cpp
//  calc
//
//  Created by Younies Mahmoud on 4/2/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

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
using namespace std;

#include "helpers.hpp"
#include "enumKTree.hpp"




Distance calculate_KN(int k , string BWT_path1 , string reverse_BWT_path1 , string BWT_path2 , string reverse_BWT_path2);

int main(int argc, const char * argv[]) {
    // insert code here...
    
    vector<string> names;
    names.push_back( "xenoRefMrna");
    names.push_back( "hu12");
    names.push_back( "Ms_lemurextended");
    
    
    string folder = "/Users/youniesmahmoud/Documents/spring2016/biosequence_algorithms/project/bio_project/";
    
    int k = 21;
    
    
    string output_folder =  "/Users/youniesmahmoud/Documents/spring2016/biosequence_algorithms/cse584a_sp16/final/third_stage/data_large/output_matrix223.csv";
    string output_distance = "/Users/youniesmahmoud/Documents/spring2016/biosequence_algorithms/cse584a_sp16/final/third_stage/data_large/output_distance_matrix223.csv";
    string output_shared = "/Users/youniesmahmoud/Documents/spring2016/biosequence_algorithms/cse584a_sp16/final/third_stage/data_large/output_shared_matrix223.csv";
    
    ofstream  puttingResult(output_folder, std::ofstream::out | std::ofstream::app);
    ofstream  puttingDistance(output_distance, std::ofstream::out | std::ofstream::app);
    ofstream  puttingShared(output_shared, std::ofstream::out | std::ofstream::app);
    
    
    
    puttingDistance << "    ,";
    puttingResult  <<"  ,";
    
    puttingShared <<  "  ,";
    
    for(int i = 0 ; i < 2 ; ++i)
    {
        puttingDistance << names[i] <<",";
        puttingResult << names[i] <<",";
        puttingResult << names[i] <<",";
        puttingResult << names[i] <<",";
        puttingShared << names[i] << ",";
    }
    
    puttingDistance << endl;
    puttingResult  << endl;
    
    puttingShared << endl;
    
    for(long long i = 0 ; i < 1 ;++i)
    {
        
        puttingDistance << names[i] <<",";
        puttingResult << names[i] <<",";
        puttingShared << names[i] << ",";
        
        for(long long j = 1 ; j < 2 ; j++)
        {
            
            string forward_1 = folder + names[i] + ".bwt";
            string reverse_1 = folder + names[i] + ".bwt2";
            
            string forward_2 = folder + names[j] + ".bwt.txt";
            string reverse_2 = folder + names[j] + ".bwt2";
            cerr << forward_2 << endl;
            
            
            
            Distance distance = calculate_KN( k , forward_1, reverse_1 ,forward_2 ,reverse_2);
            
            
            cout << distance.count << endl;
            cout << distance.first_DNA_count << endl;
            cout << distance.second_DNA_count << endl;
            
            
            double minimum = min(distance.first_DNA_count,distance.second_DNA_count);
            
            double final = -1.0/21.0 * (log(distance.count/minimum));
            
            cout << final << endl;
            
            
            if(puttingResult.is_open())
                
            {
                puttingResult <<  distance.count << "," << distance.first_DNA_count << "," << distance.second_DNA_count << "," << final << ","  ;
                puttingDistance   << final << "," ;
                puttingShared << distance.count  << ",";
            }
        }
        
        
        puttingDistance << endl;
        puttingResult << endl;
        puttingShared << endl;
        
    }
    puttingDistance.close();
    puttingResult.close();
    puttingShared.close();
    return 0;
}




Distance calculate_KN(int k , string BWT_path1 , string reverse_BWT_path1 , string BWT_path2 , string reverse_BWT_path2)
{
    BWT2 first_bwt;
    BWT2 second_bwt;
    
    
    ifstream BWT_file(BWT_path1);
    ifstream reverse_BWT_file(reverse_BWT_path1);
    
    if (BWT_file.is_open())
    {
        getline(BWT_file,  first_bwt.BWT);
        
        BWT_file.close();
    }
    else
    {
        cerr << "no file " << endl;
    }
    
    if (reverse_BWT_file.is_open())
    {
        getline(reverse_BWT_file,  first_bwt.BWT_dag);
        
        reverse_BWT_file.close();
    }
    
    
    //first_bwt.BWT = "AAAAAAAAA$";
    //first_bwt.BWT_dag = "AAAAAAAAA$";
    
    
    
    ifstream BWT_file2(BWT_path2);
    ifstream reverse_BWT_file2(reverse_BWT_path2);
    
    if (BWT_file2.is_open())
    {
        getline(BWT_file2,  second_bwt.BWT);
        
        BWT_file2.close();
    }
    else
    {
        cerr << "no file " << endl;
    }
    
    if (reverse_BWT_file2.is_open())
    {
        getline(reverse_BWT_file2,  second_bwt.BWT_dag);
        
        reverse_BWT_file2.close();
    }
    
    
    //second_bwt.BWT = "AAAAAAAAA$";
    //second_bwt.BWT_dag = "AAAAAAAAA$";
    
    buid_bwt(first_bwt);
    buid_bwt(second_bwt);
    
    
    
    
    Distance distance = enum_difference( k , first_bwt , second_bwt );
    
    cout << distance.count<< endl;
    cout << distance.first_DNA_count << endl;
    cout << distance.second_DNA_count << endl;
    
    
    return distance;
    
}




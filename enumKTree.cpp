//
//  enumKTree.cpp
//  calculating_distance
//
//  Created by Younies Mahmoud on 5/2/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "enumKTree.hpp"






Distance enum_difference(int k , BWT2 &DNA1 , BWT2 &DNA2 )
{
    Distance all_counts;
    all_counts.count = 0;
    all_counts.first_DNA_count = 0;
    all_counts.second_DNA_count = 0;
    
    enumNode first_node_DNA1;
    enumNode first_node_DNA2;
    
    first_node_DNA1.I.first = 0;
    first_node_DNA1.I.second = (int)DNA1.BWT.size() -1;
    first_node_DNA1.I_dag.first = 0;
    first_node_DNA1.I_dag.second = (int)DNA1.BWT.size() - 1;
    first_node_DNA1.depth = 0;
    
    first_node_DNA2.I.first = 0;
    first_node_DNA2.I.second = (int)DNA2.BWT.size() - 1;
    first_node_DNA2.I_dag.first = 0;
    first_node_DNA2.I_dag.second = (int)DNA2.BWT.size() - 1;
    first_node_DNA2.depth = 0;
    
    
    stack< enumNode > stack_first_DNA;
    stack_first_DNA.push(first_node_DNA1);
    
    stack< enumNode > stack_second_DNA;
    stack_second_DNA.push(first_node_DNA2);
    
    while (!stack_first_DNA.empty())
    {
 //    cerr << "dodod!\n";
        if(stack_second_DNA.empty())
            cerr << "verrry verry horrible error !!! \n";
        
        enumNode first_DNA_node = stack_first_DNA.top();
        stack_first_DNA.pop();
        enumNode second_DNA_node = stack_second_DNA.top();
        stack_second_DNA.pop();
        
        if (first_DNA_node.depth == k)
        {
 //           cerr << "here2\n";
            if(has_interval(first_DNA_node))
                all_counts.first_DNA_count += 1;
            if(has_interval(second_DNA_node))
                all_counts.second_DNA_count += 1;
            if(has_interval(first_DNA_node) && has_interval(second_DNA_node))
                all_counts.count += 1;
        }
        else
        {
            for (int base = 1 ; base < 5 ; ++base)
            {
                enumNode first_child;
                enumNode second_child;
                char base_char = base + '0';
                
                bool first_flag = computing_child(first_DNA_node, base_char, DNA1, first_child);
                
                bool second_flag = computing_child(second_DNA_node, base_char, DNA2, second_child);
                
                
                if((!first_flag) && (!second_flag))
                {
                    //cerr << "hoohho 2\n";
                    continue;
                }
                
                stack_first_DNA.push(first_child);
                stack_second_DNA.push(second_child);
            }
        }
        
    }
    
    return all_counts;
    
}



bool computing_child(enumNode &parent, char base , BWT2 &DNA , enumNode &newChild )
{
    newChild.depth = 0;
    newChild.depth = parent.depth + 1;
    if(!has_interval(parent))
    {
        newChild.I = parent.I;
        newChild.I_dag = parent.I_dag;
        return  false;
    }
    
    pair<int, int> newI;
    
    if(parent.I.first != 0)
    {
        newI.first  = DNA.C[base - '0'] + DNA.occ_B[base - '0'][ parent.I.first - 1];// + 1;//i ← C[P[k]] + occB(P[k], i − 1) + 1
        newI.second = DNA.C[base - '0'] + DNA.occ_B[base - '0'][parent.I.second] - 1;//j ← C[P[k]] + occB(P[k], j)
    }
    else
    {
        newI.first  = DNA.C[base - '0'] ;//+ 1;//i ← C[P[k]] + occB(P[k], i − 1) + 1
       // cerr << base - '0' << endl;
        //cerr <<DNA.occ_B[base - '0'][parent.I.second] - 1 << endl;
        newI.second = DNA.C[base - '0'] + DNA.occ_B[base - '0'][parent.I.second] - 1;//j ← C[P[k]] + occB(P[k], j)
    }
      // if(newI.first >= newI.second )
        //return false;
    
    
    // computing i' , j'
    pair<int, int> newI_dag;
    
    //to calculate x
    int x = 0;
    for (int i = 0 ; i < base - '0' ; ++i)
    {
        if(parent.I.first == 0)
            x += DNA.occ_B[i][parent.I.second];
        else
            x += DNA.occ_B[i][parent.I.second] - DNA.occ_B[i][parent.I.first - 1];
        
    }
    
    //to calculate y
    int y;
    if(parent.I.first != 0)
        y = x + DNA.occ_B[base - '0'][parent.I.second] - DNA.occ_B[base - '0'][parent.I.first - 1];
    else
        y = x +  DNA.occ_B[base - '0'][parent.I.second];
    
    
    
    newI_dag.first = parent.I_dag.first + x;
    newI_dag.second = parent.I_dag.first + y - 1;
    
    
    newChild.I = newI;
    newChild.I_dag = newI_dag;
    newChild.depth = parent.depth + 1;
    
  
    

    return has_interval(newChild);
    
}







bool has_interval(enumNode node){
    return node.I.second >= node.I.first;
}












void buid_bwt(BWT2 & DNA)
{
 
    
    
    int bwt_length = (int)DNA.BWT.size();
    
    vector<int> C(5,0);
    DNA.C = C;
    
    
    for (int i = 0 , n = bwt_length ; i < n ; ++i)/// building c arrays  and for making A C G T --> 1 2 3 4
    {
        DNA.BWT[i]          =   convert_DNA_to_number(DNA.BWT[i]);
       // DNA.BWT_dag[i]  =   convert_DNA_to_number(DNA.BWT_dag[i]);
        if(DNA.BWT[i] - '0' + 1 >= 5)
            continue;
        DNA.C[DNA.BWT[i] - '0' + 1 ]++;
    }
    
    for (int i = 1 , n = (int)C.size();  i < n  ; ++i)// to complete C's
        DNA.C[i] += DNA.C[i -1];
    
    
    
    
    
    
    DNA.occ_B =  vector<vector<int> > ( 5 , vector<int>(bwt_length , 0));// for occ in BWT
//DNA.occ_reverse_B =  vector<vector<int> > ( 5 , vector<int>(bwt_length , 0));// for occ in reverse BWT
    
    for (int i = 0 , n = (int)DNA.occ_B.size(); i < n; ++i)
    {
        for(int j = 0 ; j < bwt_length ; ++j)
        {
            if(DNA.BWT[j] - '0' == i)
            {
                if(j > 0)
                    DNA.occ_B[i][j] = DNA.occ_B[i][j - 1] + 1;
                else
                    DNA.occ_B[i][j]++;
            }
            else if(j > 0)
                DNA.occ_B[i][j] = DNA.occ_B[i][j - 1];
            
        /*    if (DNA.BWT_dag[j] - '0' == i)
            {
                if(j > 0)
                    DNA.occ_reverse_B[i][j] = DNA.occ_reverse_B[i][j - 1] + 1;
                else
                    DNA.occ_reverse_B[i][j]++;
            }
            else if(j > 0)
                DNA.occ_reverse_B[i][j] = DNA.occ_reverse_B[i][j - 1];
         */
            
        }
    }
    
    
    
    
}





















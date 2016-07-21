//
//  helpers.cpp
//  calculating_distance
//
//  Created by Younies Mahmoud on 5/2/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//


#include "helpers.hpp"






char convert_DNA_to_number(char base)
{
    switch (base) {
        case 'A':
            return '1';
            break;
        case 'C':
            return '2';
            break;
        case 'G':
            return '3';
            break;
        case 'T':
            return '4';
            break;
        case '$':
            return '0';
            break;
            
        default:
            cerr << "error" << endl;
    }
    
    return 'E';
}
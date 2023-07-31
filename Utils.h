#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <stdbool.h>
#include <bitset>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <sstream>

//outil pour print des msg de debug
#define DEBUG(msg) cout<<"[DEBUG] "<<msg<<endl

using namespace std;

/**
 * Structure allowing to connect 2 residue. Used in the reading of the BLOSUM matrix
*/
struct BlosumCell{

    char residue1;
    char residue2;

    bool operator==(const BlosumCell &o) const {
        return (o.residue1 == residue1 && o.residue2 == residue2);
    }

    bool operator<(const BlosumCell &o) const{
        return (o.residue1 < residue1 || (o.residue1 == residue1 && o.residue2 < residue2));
    }
};


/**
 * Function that swap from litle to big endian
*/
int32_t swap_int32(int32_t val);


/**
 * Function that receive the number of the residue and @return the residue associated
*/
char decodeResidue(int value);


/**
 * Function that reads the file and get a list of offset (result)
*/
void get_offset_table(ifstream* file,int32_t length, int32_t* result);


/**
 * Function that reads the BLOSUM file (path) and create a matrix map easyer to use for the client
*/
void init_blosum_matrix(string path, map<BlosumCell,int>* map);

#endif
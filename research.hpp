#ifndef RESEARCH_HPP
#define RESEARCH_HPP

#include <string>
#include <fstream>
#include <iostream>
#include <stdbool.h>
#include <bitset>
#include "Utils.h"
#include "Database.h"
#include "Tas_min.hpp"



typedef vector<vector<int>> Matrix;

/**
 * Function that calculate the score of each protein with the Rognez amelioration
*/
int Smith_Waterman(Database* db,string seq1, string seq2,map<BlosumCell,int>* blosum);


/**
 * Smith_waterman, longer and heavier in term of memory
*/
int Smith_Waterman2(Database* db,string seq1, string seq2, map<BlosumCell,int>* blosum);

/**
 * Function that open the fasta.psq file and reads it
**/
bool search_top_20(Database* db, Protein* requestprot, Tas_min* top_20);

/**
 * Function that reads the fasta file and add the sequence to the prot
*/
void get_protein_from_request_file(string file, Protein* prot);


/**
 * Get the headers of all the Proteins of the list and prints their header and score
*/
void get_headers(Protein proteins[20], Database* db);


/**
 * Function that prints the header of the protein
*/
void print_header(Protein* prot, string result);


/**
 * Function that transform the tas_min to a list of protein with the first element that is the min score
 * If the score are the same, sort by the iterator of offset index (from highest to lowest)
 * Modify list of Protein received by reference sorted by score (from lowest to highest)
*/
void make_increasing_list(Tas_min* top20, Protein prot[20]);


/**
 * Projet_prelim
 * Function that search the same protein and set the offset number  in index of the request_prot
*/
bool search_protein(Database* db, Protein* request_prot);

/**
 * Projet_prelim
 * Function that gets the header of the protein request_prot
*/
void get_header(Database* db, Protein* request_prot);

/**
 * Projet_prelim
 * Function that compares 2 sequences 
 * @return true if they are the same
*/
bool compare_sequences(string seq1, string seq2);



/**
 * Function that calculate the score of the seq
 * Add the protein in the top_20 if score > min
 * @return the new minimum value
*/
int calculate_score(int min, string seq, Tas_min* top_20, Database* db, Protein* request_prot, map<BlosumCell, int>*  blosum, int i);
#endif
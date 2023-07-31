#ifndef RESEARCHOPT_HPP
#define RESEARCHOPT_HPP

#include <thread>
#include <mutex>
#include "research.hpp"

struct thread_val{
    Database* db;
    int32_t start;
    int32_t size;
    Protein* requestprot;
    Tas_min* top_20;
    map<BlosumCell,int>* blosum;
    ifstream* fichier;
    int* min;
};

/**
 * Function that open the fasta.psq file and reads it
**/
void search_top_20opt(Database* db, Protein* requestprot, Tas_min* top_20);


/**
 * Thread handler
*/
void *read_seq(void* args);


/**
 * Function that calculate the score of the seq
 * Add the protein in the top_20 if score > min
 * uses 7 threads in parallels.
 * @return the new minimum value
*/
int calculate_score_opt(int min, string seq, Tas_min* top_20, Database* db, Protein* request_prot, map<BlosumCell, int>*  blosum, int i);
#endif
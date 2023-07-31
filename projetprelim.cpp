#include "research.hpp"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[]){
    const char * prot_path = argv[1];
    const char * db_path = argv[2];
    // Initialisation of the Database and the request protein
    Database db(db_path);
    Protein request_prot = Protein();
    get_protein_from_request_file(prot_path, &request_prot);
    search_protein(&db, &request_prot);
    get_header(&db, & request_prot);
    return 0;
}
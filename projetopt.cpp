#include "researchopt.hpp"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[]){
    char* endp = NULL;
    // Récupération des paramètres
    const char * prot_path = argv[1];
    const char * db_path = argv[2];
    const char * blosum_path = argv[3];
    int gap_open_penality = strtol(argv[4], &endp, 10);
    int gap_extension_penality = strtol(argv[5], &endp, 10);

    // Chargement de l'objet Database
    Database db(db_path);
    db.set_blosum_path((string) blosum_path);
    db.set_gap_open( gap_open_penality);
    db.set_gap_extension( gap_extension_penality);

    // Création de la protéine de recherche et du Tas min pour la recherche
    Protein request_prot = Protein();
    Tas_min top_20 = Tas_min();
    get_protein_from_request_file(prot_path, &request_prot);

    // Calcul des scores de toutes les protéines et mise à jour en continu du top20
    search_top_20opt(&db, &request_prot, &top_20);

    // Trier la liste top20
    Protein prot_list[20];
    make_increasing_list(&top_20, prot_list);

    // Récupérer et print les différents headers
    get_headers(prot_list, &db);

    return 0;
}
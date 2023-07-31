#include "research.hpp"



void get_protein_from_request_file(string file, Protein* prot){
    ifstream myfile(file);
    string line, DNA_sequence;
    if (!myfile.good()){
        cout << "Could not open file: " << file <<endl;
        return;
    }
    while (getline(myfile, line)) {
        if (line[0] == '>') {
            line.erase(0,1);
        }else {
            DNA_sequence+=line;
        }
    }
    DNA_sequence[DNA_sequence.length()] = '\0';
    prot->set_sequence(DNA_sequence);
    myfile.close();
    return;
}


bool search_top_20(Database* db, Protein* requestprot, Tas_min* top_20){
    string path = db->getPath(); //chemin de la db
    Index ind = db->getIndex();
    ifstream fichier(path+".psq", ios::binary); //fichier .psq
    int32_t* offset = ind.get_seqOffset(); //offsets de chaque séquence, stockés dans l'index
    int min = 0, intbuf = 0;
    Protein prot;
    map<BlosumCell,int> blosum;
    init_blosum_matrix(db->get_blosum_path(),&blosum); // Transforming the blosum in matrix Blosum 
    if (!fichier.good()){
		cout << "File could not be opened!" << endl;
		return false;
    }
    int32_t length = ind.get_seqNber();
    for(int32_t i = 0; i < length; i++){//pour chaque séquence
        string seq = "";
        int32_t seqSize = (offset[(int) i+1] - offset[(int) i]) -1; //taille de la séquence à lire
        while((int32_t) seq.length() < seqSize){ //lecture des charactères de la séquence
            fichier.read((char*) &intbuf, 1); //lecture d'un octet
            char residue = decodeResidue(intbuf); //conversion valeur -> charactère
            if(residue != 0){ //vérifie si le residu décodé n'est pas un bit null
                seq += residue;
            }
        }
        if (seq != ""){
            min = calculate_score(min, seq, top_20, db, requestprot, &blosum, i);
        }
    }
    fichier.close();
    return true;
}

void get_headers(Protein proteins[20], Database* db){
    ifstream header((db->getPath()) + ".phr", ios::binary);
    string result = "";
    char buffer;
    int offset, position_of_offset, stringLength;
    int32_t* offset_headers = db->getIndex().get_headerOffset();
    if(!header.good()){
        cout << "Impossible d'ouvrir le header !" << endl;
        return;
    }for (int i = 19; i >= 0; i--){
        result = "";
        position_of_offset = proteins[i].get_offset_number(); // Récupérer la position de l'offset dans la liste d'offset
        offset = offset_headers[position_of_offset];
        header.seekg(offset); // décalage vers le bon offset
        do{ //Lecture jusqu'à tomber sur 1A
            header.get(buffer); 
        }while(buffer != 0x1A);
        //parsing visibleString
        header.get(buffer); //Lecture du prochain octet
        bitset<8> size(buffer); //transformation en 8 bits distincts
        size.flip(0); //mise du 1er bit à 0
        stringLength = (int) size.to_ulong();
        char msg[stringLength];
        //lecture du visibleString
        header.read((char*) &msg, stringLength);
        result+=msg;
        print_header(&proteins[i], result);
    }
}

void print_header(Protein* prot, string result){
    //tronque la chaîne au 1er espace
    size_t pos = result.find(' ');
    result = result.substr(0, pos);
    cout << result;
    if (prot->get_score() == 0){
        cout << endl;
    }else{
        cout <<" " << prot->get_score()<<endl;
    }
}

void make_increasing_list(Tas_min* top20, Protein list[20]){
    for (int i=0; i < 20; i++){ // Sorting the list from lowest to highest score
        list[i] = top20->popMin();
    }
    for (int i=0; i<19;i++){ // If the score is the same, sorting from highest to lowest offset_number
        if (list[i].get_score() == list[i+1].get_score()){
            int j=i;
            while (list[j].get_score()==list[j+1].get_score()){
                j++;
            }
            for (int k = i; k<=j;k++){
                int z = k;
                while ((z > i) && (list[z].get_offset_number() > list[z-1].get_offset_number())){
                    swap(list[z],list[z-1]);
                    z--;
                }
            }
            i=j;
        }
    }
}


int calculate_score(int min, string seq, Tas_min* top_20, Database* db, Protein* request_prot, map<BlosumCell, int>*  blosum, int i){
    int score = Smith_Waterman(db, request_prot->get_sequence(), seq, blosum);
    if (score >= min || top_20->size() < 20){// Remplacer si le score de la protéine est supérieur à un des scores du top 20
        Protein prot = Protein(seq, score, i);
        if (top_20->size() >= 20){
            top_20->delMin();
        }
        top_20->insert_sort(prot);
        min = top_20->getMin().get_score();
    }
    return min;
}

int Smith_Waterman(Database* db,string seq1, string seq2, map<BlosumCell,int>* blosum){
    int score = 0;
    int m = seq1.length(),n = seq2.length();
    int R = db->get_gap_extension(), Q = db->get_gap_open() + R; //matrices à récup en paramètre
    int H[m],E[m],N[m],F;
    //Initialisation des tableaux
    BlosumCell cell;
    for(int i = 0; i < m; i++){
        H[i] = 0;
        E[i] = 0;
    }
    F = 0;
    for(int j = 0; j < n; j++){
        cell.residue1 = seq2.at(j);
        for(int i = 0; i < m; i++){
            cell.residue2 = seq1.at(i);
            //pour chaque cellule
            H[i] += blosum->at(cell);
            H[i] = max(H[i], F);
            H[i] = max(H[i], E[i]);
            H[i] = max(H[i], 0);
            score = max(score, H[i]);
            F -= R;
            E[i] -= R;  
            if (i+1<m) N[i+1] = H[i]; //mise en mémoire de H[i+1][j+1] dans N[i+1]
            H[i] -= Q;
            E[i] = max(H[i], E[i]);
            F = max(H[i], F);
        }
        copy(N, N+m, H); // copie de N dans H
    }
    return score;
}

int Smith_Waterman2(Database* db,string seq1, string seq2, map<BlosumCell,int>* blosum){
    int max_score = 0;
    int m = seq1.length(),n = seq2.length();
    int R = db->get_gap_extension(), Q = db->get_gap_open() + R;
    BlosumCell cell;
    //Initialisation des tableaux
    Matrix H(m+1, vector<int>(n+1)),E(m+1, vector<int>(n+1)),F(m+1, vector<int>(n+1));
    for(int j = 0; j <= n; j++){
        H[0][j] = 0;
        E[0][j] = 0;
    }for(int i = 0; i <= m; i++){
        H[i][0] = 0;
        F[i][0] = 0;
    }for(int j = 1; j <= n; j++){
        cell.residue1 = seq2.at(j-1);
        for(int i = 1; i <= m; i++){
            cell.residue2 = seq1.at(i-1);
            E[i][j] = max(H[i][j-1] - Q, E[i][j-1] - R);
            F[i][j] = max(H[i-1][j] - Q, F[i-1][j] - R);
            H[i][j] = max(max(H[i-1][j-1]+blosum->at(cell), 0), max(E[i][j], F[i][j]));
            max_score = max(max_score, H[i][j]); 
        }
    }
    return max_score;
}


bool search_protein(Database* db, Protein* request_prot){
    string path = db->getPath(); //chemin de la db
    Index ind = db->getIndex();
    ifstream fichier(path+".psq", ios::binary); //fichier .psq
    int32_t* offset = ind.get_seqOffset(); //offsets de chaque séquence, stockés dans l'index
    int intbuf = 0;
    if (!fichier.good()){
		cout << "File could not be opened!" << endl;
		return false;
    }
    int32_t length = ind.get_seqNber();
    for(int32_t i = 0; i < length; i++){ //pour chaque séquence
        string seq = "";
        int32_t seqSize = (offset[(int) i+1] - offset[(int) i]) -1; //taille de la séquence à lire
        while((int32_t) seq.length() < seqSize){ //lecture des charactères de la séquence
            fichier.read((char*) &intbuf, 1); //lecture d'un octet
            char residue = decodeResidue(intbuf); //conversion valeur -> charactère
            if(residue != 0){ //vérifie si le residu décodé n'est pas un bit null
                seq += residue;
            }
        }
        if (seq != ""){
            if (compare_sequences(seq, request_prot->get_sequence())){
                request_prot->set_offset_number(i);
                break;
            }
        }
    }
    fichier.close();
    return true;
}


bool compare_sequences(string seq1, string seq2){
    return seq1 == seq2;
}

void get_header(Database* db, Protein* request_prot){
    ifstream header((db->getPath()) + ".phr", ios::binary);
    string result = "";
    char buffer;
    int stringLength;
    int32_t* offset_headers = db->getIndex().get_headerOffset();
    if(!header.good()){
        cout << "Impossible d'ouvrir le header !" << endl;
        return;
    }
    header.seekg(offset_headers[request_prot->get_offset_number()]); // décalage vers le bon offset
    do{ //Lecture jusqu'à tomber sur 1A
        header.get(buffer); 
    }while(buffer != 0x1A);
    //parsing visibleString
    header.get(buffer); //Lecture du prochain octet
    bitset<8> size(buffer); //transformation en 8 bits distincts
    size.flip(0); //mise du 1er bit à 0
    stringLength = (int) size.to_ulong();
    char msg[stringLength];
    header.read((char*) &msg, stringLength); //lecture du visibleString
    result+=msg;
    print_header(request_prot, result);
}
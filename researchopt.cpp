#include "researchopt.hpp"
#include <unistd.h>

pthread_mutex_t tas_min = PTHREAD_MUTEX_INITIALIZER ,seq_read = PTHREAD_MUTEX_INITIALIZER;
int32_t* list_of_offset = 0;
int thread_working = 0;

void search_top_20opt(Database* db, Protein* requestprot, Tas_min* top_20){
    string path = db->getPath(); //chemin de la db
    Index ind = db->getIndex();
    int32_t length = ind.get_seqNber();
    pthread_t tid[7]; //Création des threads
    ifstream fichier(path+".psq", ios::binary); //fichier .psq
    if (!fichier.good()){
		cout << "File could not be opened!" << endl;
		return;
    }
    //Lancement dans thread
    list_of_offset = ind.get_seqOffset(); //offsets de chaque séquence, stockés dans l'index
    Protein prot;
    map<BlosumCell,int> blosum;
    int min = 0;
    init_blosum_matrix(db->get_blosum_path(),&blosum); // Transforming the blosum in matrix Blosum 
    int32_t start = 0;
    int32_t len_of_thread = (length -(length % 7))/7;
    for(pthread_t thread : tid){
        struct thread_val* args = (thread_val *) malloc(sizeof(struct thread_val));
        args->db = db;
        args->start = start;
        args->size = length/7;
        args->requestprot = requestprot;
        args->top_20 = top_20;
        args->blosum = &blosum;
        args->fichier = &fichier;
        args->min = &min;
        if(start == 0){
            args->size  = len_of_thread + length % 7 -1;
            start += length % 7;
        }
        start += len_of_thread;
        pthread_create(&thread, NULL,read_seq,(void*) args);
    }
    while(thread_working > 0){
        sleep(5);
    }
    fichier.close();
}

int calculate_score_opt(int min, string seq, Tas_min* top_20, Database* db, Protein* request_prot, map<BlosumCell, int>*  blosum, int i){
    int score = Smith_Waterman(db, request_prot->get_sequence(), seq, blosum);
    pthread_mutex_lock(&tas_min);
    if (score >= min || top_20->size() < 20){// Remplacer si le score de la protéine est supérieur à un des scores du top 20
        Protein prot = Protein(seq, score, i);
        if (top_20->size() >= 20){
            top_20->delMin();
        }
        top_20->insert_sort(prot);
        min = top_20->getMin().get_score();
    }
    pthread_mutex_unlock(&tas_min);
    return min;
}

void *read_seq(void* args){
    int intbuf = 0;
    struct thread_val* values = (thread_val*) args;
    thread_working++;
    for(int32_t i = values->start; i <  (values->start + values->size); i++){//pour chaque séquence
        pthread_mutex_lock(&seq_read);
        string seq = "";
        values->fichier->seekg(list_of_offset[i]);
        int32_t seqSize = (list_of_offset[(int) i+1] - list_of_offset[(int) i]) -1; //taille de la séquence à lire
        while((int32_t) seq.length() < seqSize){ //lecture des charactères de la séquence
            values->fichier->read((char*) &intbuf, 1); //lecture d'un octet
            char residue = decodeResidue(intbuf); //conversion valeur -> charactère
            if(residue != 0){ //vérifie si le residu décodé n'est pas un bit null
                seq += residue;
            }
        }
        pthread_mutex_unlock(&seq_read);
        if (seq != ""){
            *(values->min) = calculate_score_opt(*(values->min), seq, values->top_20, values->db, values->requestprot, (values->blosum), i);
        }
    }
    free(args);
    thread_working--;
    return NULL;
}




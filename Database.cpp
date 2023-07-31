#include "Database.h"

Database::Database(string path) : index(path+".pin"),db_path(path){

}

Database::~Database(){
	delete index.get_headerOffset();
    delete index.get_seqOffset();
}


void Database::set_gap_open(int gap_open) {
	gap_open_penality = gap_open;
}

void Database::set_gap_extension(int gap_extension) {
	gap_extension_penality = gap_extension;
}

void Database::set_blosum_path(string path) {
	blosum_path = path;
}


Index Database::getIndex() const{
	return index;
}

string Database::getPath() const{
	return db_path;
}


int Database::get_gap_open() const{
	return gap_open_penality;
}


int Database::get_gap_extension() const{
	return gap_extension_penality;
}

string Database::get_blosum_path() const{
	return blosum_path;
}
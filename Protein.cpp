#include "Protein.h"

// Protein
Protein::Protein(string seq, int score, int32_t offset) : sequence(seq),storedScore(score),offset_number(offset){
};

Protein::~Protein(){
};

string Protein::get_header()	const{
	return header;
};

string Protein::get_sequence()	const{
	return sequence;
};


int Protein::get_score() const{
	return storedScore;
}

int32_t Protein::get_offset_number() const{
	return offset_number;
}

void Protein::set_score(int score){
	storedScore = score;
}

void Protein::set_offset_number(int number){
	offset_number = number;
}

void Protein::set_sequence(string seq){
	sequence = seq;
};

void Protein::set_header(string head){
	header = head;
};
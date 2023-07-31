#ifndef PROTEIN_HPP
#define PROTEIN_HPP

#include <vector>
#include <string>

using namespace std;

class Protein { // Protéine en général
	private:
		string header = ""; 
		string sequence = "";
		int storedScore = 0;
		int32_t offset_number = 0;
    
	public:
		Protein(string seq = "", int score = 0, int32_t offset = 0);
		~Protein();
		int get_score() const;
		string get_sequence() const;
		string get_header() const;
		int32_t get_offset_number() const;
		void set_score(int score);
		void set_offset_number(int number);
		void set_sequence(string seq);
		void set_header(string head);
		
};

#endif
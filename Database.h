#ifndef DATABASE_HPP
#define DATABAE_HPP

#include "Protein.h"
#include "Index.h"
#include <map>
#include "Utils.h"

class Database{ // liste des protéines

	private:
		Index index;
		string db_path = "database/uniprot_sprot.fasta";
		string blosum_path = "blosum/BLOSUM62";
		int gap_open_penality = 11;
		int gap_extension_penality = 1;
	public:
		Database(string path);
		~Database();
		void set_gap_open(int gap_open);
		void set_gap_extension(int gap_extension);
		void set_blosum_path(string path);
		Index getIndex() const;
		string getPath() const;
		int get_gap_open() const;
		int get_gap_extension() const;
		string get_blosum_path() const;
};

#endif
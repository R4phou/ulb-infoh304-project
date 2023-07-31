#include "Utils.h"

//int32 endian convertisser
int32_t swap_int32(int32_t val){
    val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF ); 
    return (val << 16) | ((val >> 16) & 0xFFFF);
}

//permet la convertion d'un résidu de sa forme binaire vers un charactère
char decodeResidue(int value){
    //tableau repris de la documentation sur le BLAST (partie .psq)
	char residueTable[28] =
		{ 0,'A','B','C','D','E','F',
		'G','H','I','K','L','M','N',
		'P','Q','R','S','T','V','W',
		'X','Y','Z','U','*','O','J'};
	return residueTable[value];
}

void get_offset_table(ifstream* file, int32_t length, int32_t* result){
	int32_t buffer;
	for (int32_t i = 0; i < length;i++){
		file->read((char*)&buffer,sizeof(buffer));
		buffer=swap_int32(buffer);
		result[i] = buffer;
	}
}

void init_blosum_matrix(string path, map<BlosumCell,int>* map){
	ifstream blosum_file(path);
	if(!blosum_file.good()){
		cout <<"File not good !"<<endl;
		return;
	}
	string line;
	getline(blosum_file, line);
	bool is_first_line = true;
	vector<char> buffer;
	//lecture de la matrice blosum
	while(getline(blosum_file,line)){
		if(line.at(0) != '#' && line.length() > 0){
			//ignore les lignes commençant par #
			stringstream ss;
			ss.str(line);
			char c = ' ';
			if(is_first_line){
				//c'est la ligne contenant tous les 2nd résidus
				//lecture jusqu'en fin de ligne
				do{
					ss >> c;
					if(c != ' ' && c != '\n' && c != '\0'){
						buffer.push_back(c);
					}
				}while(ss.good());
				is_first_line = false;
			}else{
				ss >> c; //récupération du 1er résidu
				for(long unsigned int i = 0; i < buffer.size()-1; i++){
					int value;
					ss >> value;
					//création de la cellule de la matrice
					BlosumCell cell;
					cell.residue1 = c;
					cell.residue2 = buffer.at(i);
					pair<BlosumCell,int> p(cell,value);
					map->insert(p);
				}
			}
		}		
	}
	blosum_file.close();
}
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <vector>

using namespace std;

#include "Index.h"

Index::Index(string path){
    ifstream myfile(path,ios::binary);
    if (myfile.is_open()){
        int32_t buffint32, bufftime32;
        int64_t buffint64;
        myfile.read((char*)&buffint32, sizeof(buffint32));
        buffint32 = swap_int32(buffint32);
        version=buffint32; //Version
        myfile.read((char*)&buffint32, sizeof(buffint32));
        buffint32 = swap_int32(buffint32);
        dbType=buffint32; //Database type
        myfile.read((char*)&buffint32, sizeof(buffint32));
        buffint32 = swap_int32(buffint32);
        titleLength=buffint32; // Title Length

        char buffchar1[titleLength];
        myfile.read((char*)&buffchar1,sizeof(buffchar1));
        buffchar1[titleLength]='\0';
        title = buffchar1; //Title
        
        myfile.read((char*)&bufftime32, sizeof(bufftime32));
        bufftime32 = swap_int32(bufftime32);
        timeLength=bufftime32; //Time stamp length

        char buffchar2[timeLength];
        myfile.read((char*)&buffchar2,sizeof(buffchar2));
        time = buffchar2;//Time Stamp
        
        myfile.read((char*)&buffint32, sizeof(buffint32));
        buffint32 = swap_int32(buffint32);
        seqNber=buffint32; //Number of sequence (N)

        myfile.read((char*)&buffint64, sizeof(buffint64));//Pas d'inversion vers big endian car il est déjà en little endian
        residueCount=buffint64; //Residue Count
    
        myfile.read((char*)&buffint32, sizeof(buffint32));
        buffint32 = swap_int32(buffint32);
        maxSeq=buffint32; //Sequence max
        headerOffset = new int32_t [seqNber+1];
        seqOffset = new int32_t [seqNber+1];
        get_offset_table(&myfile,seqNber+1 ,headerOffset);
        get_offset_table(&myfile,seqNber+1 , seqOffset);
        myfile.close();
    }
    else{
        cout << "Couldn't open file" << endl;
    }
}

Index::~Index(){
}

int32_t Index::get_version() const{
    return version;
}

int32_t Index::get_DbType() const{
    return dbType;
}

int32_t Index::get_titleLength() const{
    return titleLength;
}

string Index::get_title() const{
    return title;
}

int32_t Index::get_timelength() const{
    return timeLength;
}

string Index::get_time() const{
    return time;
}

int64_t Index::get_residueCount() const{
    return residueCount;
}

int32_t Index::get_seqNber() const{
    return seqNber;
}

int32_t Index::get_maxSeq() const{
    return maxSeq;
}

int32_t* Index::get_headerOffset() const{
    return headerOffset;
}

int32_t* Index::get_seqOffset() const{
    return seqOffset;
}


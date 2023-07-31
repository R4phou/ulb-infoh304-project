#ifndef INDEX_HPP
#define INDEX_HPP

#include <vector>
#include <cstdint>
#include <string.h>
#include <iostream>
#include "Utils.h"
using namespace std;

class Index{ // .pin
	private:
		int32_t version;
		int32_t dbType;
		int32_t titleLength; //T
		string title;
		int32_t timeLength; //S
		string time;
		int32_t seqNber; //N
		int64_t residueCount;
		int32_t maxSeq;
		// vector<int32_t> headerOffset; // offset des header dans .phr
		int32_t* headerOffset;
		int32_t* seqOffset; // offset des sequences dans .psq
	public:
		Index(string path);
		~Index();
		int32_t get_version() const;
		int32_t get_DbType() const;
		int32_t get_titleLength() const;
		string get_title() const;
		int32_t get_timelength() const;
		string get_time() const;
		int32_t get_seqNber() const;
		int64_t get_residueCount() const;
		int32_t get_maxSeq() const;
		int32_t* get_headerOffset() const;
		int32_t* get_seqOffset() const;

};

#endif
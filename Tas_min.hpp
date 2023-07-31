#ifndef TAS_MIN_HPP
#define TAS_MIN_HPP

#include <vector>
#include "Protein.h"
using namespace std;

class Tas_min{

    private:
        vector<Protein> data;

    public:
        Tas_min();
        ~Tas_min();
        bool isEmpty() const;
        int size() const;
        Protein getMin() const;
        void insert_sort(Protein prot);
        void delMin();
        Protein popMin();
        Protein popMax();
        void rearrange(int index);
};

#endif
#include "Tas_min.hpp"


Tas_min::Tas_min(){
};

Tas_min::~Tas_min(){
};

bool Tas_min::isEmpty() const{
    return data.empty();
};

int Tas_min::size()const{
    return (int) data.size();
};

Protein Tas_min::getMin() const{
    return data[0];
};

void Tas_min::insert_sort(Protein prot){
    data.push_back(prot);
    int index = (int) data.size() - 1;
    while (index > 0) {
        int parent = (index - 1) / 2;
        if (data[parent].get_score() <= data[index].get_score()) {
            break;
        }
        swap(data[parent], data[index]);
        index = parent;
    }
};

void Tas_min::delMin() {
    // Replace the root with the last element in the heap
    data[0] = data.back();
    data.pop_back();

    // Arrange all elements recursively
    rearrange(0);
};

Protein Tas_min::popMin(){
    // Save the minimum
    Protein min = data[0];
    data[0] = data.back();
    data.pop_back();

    // Arrange all elements recursively
    rearrange(0);

    return min;
};

void Tas_min::rearrange(int index){
    int left = 2 * index + 1;
    int right = 2 * index + 2;
    int min_index = index;
    if (left < (int) data.size() && data[left].get_score() < data[min_index].get_score()) {
        min_index = left;
    }
    if (right < (int) data.size() && data[right].get_score() < data[min_index].get_score()) {
        min_index = right;
    }
    if (min_index != index) {
        swap(data[index], data[min_index]);
        rearrange(min_index);
    }
};


Protein Tas_min::popMax(){ // Complexité O(n)
    int to_delete = 0;
    Protein max = Protein();
    for (int i=0; i < (int) data.size(); i++){
        if(data[i].get_score() > max.get_score()){
            max = data[i];
            to_delete = i;
        } 
    }
    data.erase(data.begin()+to_delete);
    rearrange(0);
    return max;
};
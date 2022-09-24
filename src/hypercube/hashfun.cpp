
#include "hashfun.hpp"



HashF::HashF(const uint16_t w, const uint32_t dim) 
    : t(create_rand_float_uniform(w)), dim(dim), map(new std::map<uint32_t, bool>), vArray(new std::vector<float>(dim))
{
    for(int j=0; j<dim; j++) 
        (*vArray)[j] = create_rand_float_normal();        
}

HashF::~HashF() {
    delete map;
    delete vArray;
}



bool HashF::hash(Point *pointPtr, const uint16_t w) const {
    
    double temp = internal_product<double,float>( *(pointPtr->getVector()) , *vArray );
    temp += this->t;
    temp = temp / w;
    temp = floor(temp);

    uint64_t key = static_cast<uint32_t>(temp);
    std::map<uint32_t, bool>::iterator it;

    it = map->find(key);
    if( it == map->end()) {                  //if this h value is not already generated
        bool bit = create_rand_bernoulli();             //create a random val and insert it into map
        map->insert({key, bit});
        return bit;
    }       
    else
        return it->second;
}

HashFs::HashFs(const uint16_t k, const uint16_t w, const uint32_t tableSize, const uint32_t dim) 
    :   k(k) , w(w), tableSize(tableSize), HashFunctions(new std::vector<HashF *>(k))
{
    for(int j=0; j<k; j++) 
        (*HashFunctions)[j] = new HashF(w, dim);   
}

HashFs::~HashFs() {
    for(int i=0; i<HashFunctions->size(); i++)
        delete (*HashFunctions)[i];
    delete HashFunctions;
}

uint32_t HashFs::amplified_hash(Point *pointPtr) const {

    uint32_t key = 0;
    int i,j,z;

    // calculate the F_i(h_i(p)))  hashfunctions values
    for(i=0; i< this->k; i++) {
        bool bit = (*HashFunctions)[i]->hash(pointPtr, this->w);
        uint32_t bitInt = (bit) ? static_cast<uint32_t>(1) : static_cast<uint32_t>(0);      // convert it to uint32_t
        bitInt = bitInt << i;          // shift left bit value   
        key = key | bitInt; 
        //hashFunctions_vals.push_back(bit);
    }
    /*
    uint32_t key = 0;
    uint16_t size = hashFunctions_vals.size();
    for( i=size-1; i >= 0; --i ) {
        bool hashed_val = hashFunctions_vals[i];                    // get the hashed value as boolean, with reverse order
        
    }*/

    return key;
}

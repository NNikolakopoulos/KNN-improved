#include "hashf.hpp"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   CONSTRUCTOR  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


HashFunctions::HashFunctions(const uint16_t k, const uint16_t w, const uint16_t M, const uint16_t L,const uint32_t tableSize, const uint32_t dim) 
    : k(k), w(w), M(M), L(L),tableSize(tableSize), dim(dim)
{
    int i;
    float rand_float;
    std::vector<float> *tempV, *tempT;
    std::vector<std::vector<float> *> *tempVector;

    vArray = new std::vector<std::vector<std::vector<float>*>*>;
    rArray = new std::vector<uint32_t>;
    tArray = new std::vector<std::vector<float>*>;

    for(int z=0; z<this->L; z++) {
        tempVector = new std::vector<std::vector<float> *>;
        tempT = new std::vector<float>;
        for(i=0; i<k; i++) {
            tempV = new std::vector<float>;
            for(int j=0; j<this->dim; j++) {
                tempV->push_back(create_rand_float_normal());
            }
            tempVector->push_back(tempV);                                // push v_i vector for the h_i hashfunction
            //tempV.clear();
            tempT->push_back(create_rand_float_uniform(w));              // push t_i number for the h_i hashfunction
            if(z == 0)                                                  // we need to create only one r_i vector
                rArray->push_back(create_rand_int_uniform(UINT32_MAX));
        }
        vArray->push_back(tempVector);                                   // push all the h hashfunctions of this particular g hashfunction.
        //tempVector.clear();
        tArray->push_back(tempT);
        //tempT.clear();
    }
}

HashFunctions::~HashFunctions() {
    delete rArray;
    for(int i=0; i< vArray->size(); i++) {
        delete (*tArray)[i];
        for(int j=0; j< (*vArray)[i]->size(); j++)
            delete (*(*vArray)[i])[j] ;
        delete (*vArray)[i];
    }
    delete tArray;
    delete vArray;
}


// create a random float number with normal distribution
float create_rand_float_normal() {
    std::default_random_engine generator;
    std::normal_distribution<float> distribution (0.0,1.0);
    return distribution(generator);
}

// create a random float number with uniform distribution
float create_rand_float_uniform(const uint32_t w) {
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution (0.0,w);
    return distribution(generator);
}

float create_rand_int_uniform(const uint32_t w) {
    std::default_random_engine generator;
    std::uniform_int_distribution<uint32_t> distribution (0.0,w);
    return distribution(generator);
}

bool create_rand_bernoulli() {
    std::default_random_engine generator;
    std::bernoulli_distribution distribution (0.5);
    return distribution(generator);
}

float create_rand_double_uniform(const double window) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution (0.0,window);
    return distribution(generator);
}






/* ~~~~~~~~~~~~~~~~~~~~~~   hashfunctions h  ~~~~~~~~~~~~~~~~~~~~~~~ */


uint64_t HashFunctions::amplified_hash(Point *pointPtr, const uint16_t g_i) const {

    std::vector<double> h_hashFunctions;
    int i,j,z;

    // calculate the h-hashfunctions values
    for(i=0; i< this->k; i++) {
        double temp = internal_product<double,float>( *(pointPtr->getVector()) , (*(*(*vArray)[g_i])[i]));
        temp += (*(*tArray)[g_i])[i];
        temp = temp /this->w;
        temp = floor(temp);
        h_hashFunctions.push_back(temp);
    }
    uint64_t sum = 0;
    for(i=0; i< this->k; i++) {
        uint64_t temp = (*rArray)[i] * h_hashFunctions[i];
        sum += mod(temp, this->M);
    }
    sum = mod(sum, this->M);                    // this is the hashed ID of the point

    return sum;
    //pointPtr->setHashedId(sum);

    //return mod(sum, this->tableSize);            // this is the amplified function final value
}


uint64_t HashFunctions::amplified_hash(std::vector<double> *vector, const uint16_t g_i) const {

    std::vector<double> h_hashFunctions;
    int i,j,z;

    // calculate the h-hashfunctions values
    for(i=0; i< this->k; i++) {
        double temp = internal_product<float,double>( (*(*(*vArray)[g_i])[i]), *vector , (vector->size())/2);
        temp += (*(*tArray)[g_i])[i];
        temp = temp /this->w;
        temp = floor(temp);
        h_hashFunctions.push_back(temp);
    }
    uint64_t sum = 0;
    for(i=0; i< this->k; i++) {
        uint64_t temp = (*rArray)[i] * h_hashFunctions[i];
        sum += mod(temp, this->M);
    }
    sum = mod(sum, this->M);                    // this is the hashed ID of the point

    return sum;
    //pointPtr->setHashedId(sum);

    //return mod(sum, this->tableSize);            // this is the amplified function final value
}
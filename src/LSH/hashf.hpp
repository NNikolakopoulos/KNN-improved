//#ifndef HASHF_HEADERS
//#define HASHF_HEADERS
#pragma once
#include "../../headers.hpp"
#include "hash.hpp"
#include "../operations.hpp"

float create_rand_float_normal();
float create_rand_float_uniform(const uint32_t window);
float create_rand_double_uniform(const double window);
float create_rand_int_uniform(const uint32_t window);
bool create_rand_bernoulli();

class HashFunctions {
    private:
        std::vector<std::vector<std::vector<float>*>*> *vArray;
        std::vector<uint32_t> *rArray;
        std::vector<std::vector<float>*> *tArray;
        uint16_t k;
        uint16_t w;
        uint16_t M;
        uint16_t L;
        uint32_t tableSize;
        uint32_t dim;

    public:
    HashFunctions(const uint16_t k, const uint16_t w, const uint16_t M, const uint16_t L,const uint32_t tableSize, const uint32_t dim); 
    ~HashFunctions();

    uint64_t amplified_hash(Point *pointPtr, const uint16_t g_i) const;
    uint64_t amplified_hash(std::vector<double> *vector, const uint16_t g_i) const;
};

//#endif
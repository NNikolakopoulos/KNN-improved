//#ifndef OPERATIONS_HEADERS
//#define OPERATIONS_HEADERS
#pragma once
#include "../headers.hpp"
#include "./LSH/hash.hpp"

template<typename T,typename U> 
double internal_product(std::vector<T> &a, std::vector<U> &b) {
    double sum = 0;
    for(int i=0; i<a.size(); i++) {
        sum = sum + (a[i] * b[i]);

    }
    
    return sum;
}
template<typename T,typename U> 
double internal_product(std::vector<T> &a, std::vector<U> &b, const int numDimensions) {
    double sum = 0;
    for(int i=0; i<numDimensions; i++) {
        sum = sum + (a[i] * b[i]);

    }
    
    return sum;
}
uint32_t pow2(const uint32_t num);
uint64_t mod(const uint64_t a, const uint64_t b);
double euclideanDistance(std::vector<double> &a, std::vector<double> &b);
double euclideanDistance1D(const double a,const double b);
uint32_t hammingDistance(const uint64_t n1, const uint64_t n2);
std::queue<uint32_t> & hamDist1(std::queue<uint32_t> & queue,  const uint32_t key, std::unordered_set<uint32_t> & expl_set_hamDist, const uint32_t length);
int search_num_range(std::vector<double> &vec, double &num);
double compute_delta();
double discrete_frechet(Curve *c1, Curve *c2);
double optimal_traversal(Curve *c1, Curve *c2);


//#endif
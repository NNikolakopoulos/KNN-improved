#pragma once
#include "../../headers.hpp"
#include "hyperc.hpp"
#include "../operations.hpp"


class HashF {
    private:
        float t;
        uint32_t dim;
        std::map<uint32_t, bool> *map;
        std::vector<float> *vArray;
        
    public:
        HashF(const uint16_t w, const uint32_t dim); 
        ~HashF();

        bool hash(Point *pointPtr, const uint16_t w) const;
};

class HashFs {
    private:
        uint16_t k;
        uint16_t w;
        uint32_t tableSize;
        std::vector<HashF *> *HashFunctions;
        
    public:

        HashFs(const uint16_t k, const uint16_t w, const uint32_t tableSize, const uint32_t dim); 
        ~HashFs();

        uint32_t amplified_hash(Point *pointPtr) const;
};

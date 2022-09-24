#pragma once

#include "../../headers.hpp"
#include "../operations.hpp"
#include "../inputHandle.hpp"
#include "../LSH/hash.hpp"
#include "hashfun.hpp"


/*  Point is declared in LSH/hash.hpp */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HYPERCUBE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

class HyperCube {
    private:
        std::vector<std::list<Point *> *> *array;
        std::unordered_map<std::string, uint64_t> ids;
        uint64_t size;

    public:
        HyperCube(const uint32_t size);
        ~HyperCube();

        const int insert(Point *pointPtr, const uint64_t key);

        std::vector<std::list<Point *> *> *getArray() { return this->array; }

        std::multimap<double, std::pair<Point *, std::chrono::duration<double>> > & bruteForce(Point *pointPtr, const uint32_t K,std::multimap<double, std::pair<Point *, std::chrono::duration<double>> > & map) const;
        void KNN(Point *pointPtr, const u_int32_t key, const uint32_t K,std::multimap<double,std::pair<Point *, std::chrono::duration<double>>> & map,  userInput *ui) const;
        int rangedSearchAmplified(Point *pointPtr, const uint64_t key, std::unordered_set<Point *> *point_set, userInput *ui, const double R, const bool addBarrier) ;
        
};

class RandomizedProjection {
    private:
        userInput *ui;
        HashFs *hashFun;
        HyperCube *hyperCube;
        

    public:
        RandomizedProjection(std::vector<Point *> * PointsVector, userInput *ui);
        ~RandomizedProjection();

        std::multimap<double, std::pair<Point *, std::chrono::duration<double>>> queryKNN(Point *pointPtr, const uint32_t K) const;

        bool queryFile(std::string queryFile) const;

        int rangedSearchCentroid(Point *point, std::unordered_set<Point *> *cluster, const double R, const bool addBarrier);
        void unassignedPoints(std::unordered_map<Point *, std::unordered_set<Point *> *> *clusters, std::unordered_set<Point *> *centroids);


};
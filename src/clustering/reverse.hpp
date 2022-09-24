#pragma once

#include "../../headers.hpp"
#include "../inputHandle.hpp"
#include "../LSH/hash.hpp"
#include "../hypercube/hyperc.hpp"

class reverseAssignment {
    private:
        LSH *lsh;
        RandomizedProjection *cube;
        LSH_discrete_frechet *lsh_frechet;
        std::unordered_set<Point *> *centroids;
        int method;                                    // defined in headers.hpp

    public:
        reverseAssignment(std::vector<Point *> *Points, userInput *ui, std::unordered_set<Point *> *centroids, const int method);       
        ~reverseAssignment();

        double calculateStartRadius() const;
        std::unordered_map<Point *, std::unordered_set<Point *> *> *assignment(std::unordered_map<Point *, std::unordered_set<Point *> *> *clusters, 
            std::unordered_set< Point *> *centroids);

};
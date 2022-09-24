/*#ifndef HASH_HEADERS
#define HASH_HEADERS*/
#pragma once

#include "../../headers.hpp"
#include "hashf.hpp"
#include "../operations.hpp"
#include "../inputHandle.hpp"

#include "../fred/types.hpp"
#include "../fred/curve.hpp"
#include "../fred/frechet.hpp" 

fred::Curve * convertCurve(std::vector<double> *filtered);
void assignPoint(std::unordered_map<Point *, std::unordered_set<Point *> *> *clusters, Point *point, const int typeFrechet);
void assignPointCloserCluster(std::unordered_map<Point *, std::unordered_set<Point *> *> *clusters, std::list<Point *> & concurrentCentroids, Point *point, const int typeFrechet);


class Point {
    private:
        std::vector<double> *vect;
        std::string itemId;

    public:
        Point(std::string id, std::vector<double> *vector);
        Point(std::tuple<std::string, std::vector<double> *> tuple);
        ~Point();
        
        std::vector<double> *getVector() const { return vect; }
        std::string getItemId() const { return itemId; }

        inline void add(const double i) { this->vect->push_back(i); }
        inline void setVector(std::vector<double> *newVect) { this->vect = newVect; }

        inline void deleteVector() { delete vect; }

};


class HashTable {
    private:
        std::vector<std::list<Point *> *> *array;
        std::unordered_map<std::string, uint64_t> ids;
        uint64_t size;
 
    public:
        HashTable(const uint32_t tableSize);
        ~HashTable();

        const int insert(Point *pointPtr, const uint64_t key);
        void deletePoints();

        std::vector<std::list<Point *> *> * getArray() { return this->array; }

        std::tuple<Point *, double> NN(Point *pointPtr, const uint64_t key) const;
        std::multimap<double, std::pair<Point *, std::chrono::duration<double>> > & bruteForce(Point *pointPtr, const uint32_t K,
            std::multimap<double, std::pair<Point *, std::chrono::duration<double>> > & map, int typeFrechet) const;
        void KNN(Point *pointPtr, const u_int64_t key, const uint32_t K,std::multimap<double,std::pair<Point *, std::chrono::duration<double>>> & map, 
            std::unordered_set<std::string> & explored_set, int typeFrechet) const;
        int rangedSearchAmplified(Point *pointPtr, const uint64_t key, const double R, std::unordered_set<Point *> * point_set, 
            std::unordered_set<std::string> & explored_set, const bool addBarrier, const int typeFrechet);


        
};

class LSH {
    private:
        int typeFrechet;
        userInput *ui;
        HashFunctions *hashFun;
        std::vector<HashTable *> HTarray;

    public:
    LSH(std::vector<Point *> * PointsVector, userInput *ui, int typeFrechet);
    LSH(std::vector<Point *> * PointsVector, userInput *ui, const double delta, std::pair<double,double> *t ,int typeFrechet);
    LSH(std::vector<Point *> * PointsVector, userInput *ui, const double delta, double t , double epsilon,int typeFrechet);
    ~LSH();

    const std::vector<HashTable *>& getHTarray() const { return HTarray; }
    std::multimap<double, std::pair<Point *, std::chrono::duration<double>>> &queryKNN(Point *pointPtr, const uint32_t K, // 
      std::multimap<double, std::pair<Point *, std::chrono::duration<double>> >&map, std::unordered_set<std::string> &explored_set, const double delta, std::pair<double,double> *t, const double epsilon) const;

    int rangedSearchCentroid(Point *point, std::unordered_set<Point *> *cluster, const double R, const bool addBarrier);
    int rangedSearchCentroid(Point * centroid, std::unordered_set<Point *> *cluster, const double R, const bool addBarrier,std::unordered_set<std::string> &explored_set, const double delta, std::pair<double,double>& t );
    void unassignedPoints(std::unordered_map<Point *, std::unordered_set<Point *> *> *clusters, std::unordered_set<Point *> *centroids);

    bool queryFile(std::string queryFile) const;
    
    std::vector<std::pair<double,double>*> * snap(Curve &curve, const double delta, std::pair<double,double> &t) const;
    std::vector<double> * getConcatCurve(std::vector<std::pair<double,double>*> *snapped_curve) const;

    std::vector<double> * filter(Curve &curve, const double epsilon) const;
    std::vector<double> * snap_continuous(std::vector<double> *filtered, const double delta, double t) const;
    std::vector<double> * minima_maxima(std::vector<double> *filtered) const;
};


// define shifted grid class
class Grid {
    private:
        double delta;  
        double epsilon;                             
        std::pair<double,double> *t;                      // random vector t in (0,delta)^dim(in our case dim=2), to shift grid
        LSH *lsh;                                   // each grid has each own LSH 

    public:
    Grid(std::vector<Curve *> * CurveVector, userInput *ui,const double delta);
    Grid(std::vector<Curve *> * CurvesVector, userInput *ui,const double delta, const double epsilon);
    ~Grid();

    LSH *getLSH() const {  return lsh; }

    std::multimap<double, std::pair<Curve *, std::chrono::duration<double>>> & GridQueryKNN(Curve *CurvePtr, const uint32_t K, // 
      std::multimap<double, std::pair<Curve *, std::chrono::duration<double>> >&map, std::unordered_set<std::string> &explored_set);

    bool rangedSearchCentroidGrid(Point *point, std::unordered_set<Point *> *cluster, std::unordered_set<std::string> &explored_set, const double R, const bool addBarrier, const int typeFrechet);
};

class LSH_discrete_frechet {
    private:
        userInput *ui;
        std::vector<Grid *> grids;

    public:
        LSH_discrete_frechet(std::vector<Curve *> * CurveVector, userInput *ui, const double delta);
        ~LSH_discrete_frechet();

        bool queryFile(std::string queryFile) const;

        int rangedSearchCentroidFrechet(Point *point, std::unordered_set<Point *> *cluster, const double R, const bool addBarrier);
        void unassignedPoints(std::unordered_map<Point *, std::unordered_set<Point *> *> *clusters, std::unordered_set<Point *> *centroids);
};

class LSH_continuous_frechet {
    private:
        userInput *ui;
        Grid *grid;

    public:
        LSH_continuous_frechet(std::vector<Curve *> * CurveVector, userInput *ui, const double delta, const double epsilon);
        ~LSH_continuous_frechet() { delete grid; }

        bool queryFile(std::string queryFile) const;

};
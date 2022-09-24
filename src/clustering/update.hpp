#pragma once
#include "../../headers.hpp"
#include "../operations.hpp"

bool update(std::unordered_map<Point*, std::unordered_set<Point*>*>*, std::unordered_set<Point*>*);
bool update_curves(std::unordered_map<Curve*, std::unordered_set<Curve*>*>*, std::unordered_set<Curve*>*);
std::vector<double>* mean_vectors(std::unordered_set<Point*>*);
std::list<std::tuple<double, double, double, double>> opt_traversal(Curve*, Curve*);
Curve* mean_curve(Curve*, Curve*);
Curve* mean_curve(std::unordered_set<Curve*>*);


class CBT {
    private:
        std::vector<Curve *> * array;
        int leaf_level_index;
        int size;
    
    public:
        CBT(std::unordered_set<Curve *> &array, int numLeaves);
        ~CBT() { delete array; }

        void setCurve(int index, Curve *curve)  { if(index<size) array->at(index) = curve; }
        Curve *getCurve(int index)    { return ( (index<size) ? array->at(index) : nullptr ); }

        int getLeft(int index)    { return 2*index+1; }
        int getRight(int index)   { return 2*index+2; }

        bool isLeaf(int index)    { return index >= leaf_level_index; }
        bool isInternal(int index)    { return index < leaf_level_index; }
};
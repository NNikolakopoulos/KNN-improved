#pragma once
#include "../../headers.hpp"

std::unordered_set<Point*>* initialization(std::vector<Point*>*, const unsigned int);
std::unordered_set<Point*>* initialization_curves(std::vector<Point*>*, const unsigned int);
void deleteClusters(std::unordered_map<Point*, std::unordered_set<Point*>*>*);
std::vector<double> silhouette(std::unordered_map<Point*, std::unordered_set<Point*>*>*, const int);
void deletePointsVector(std::vector<Point *> *, std::string);
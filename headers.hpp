//#ifndef HEADERS_GUARD
//#define HEADERS_GUARD
#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <list>
#include <random>
#include <bitset>
#include <fstream>
#include <cmath>
#include <map>
#include <queue>
#include <iterator>
#include <chrono>
#include <algorithm>
#include <bits/stdc++.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>

#define THRESHOLD 1000
#define PADDING 10000
#define TYPE_VECTOR 0
#define TYPE_DISCRETE 1
#define TYPE_CONTINUOUS 2
#define TYPE_CUBE 3

class Point;
class HashTable;
class LSH;
class HashFunctions;
class HashFs;
typedef Point Curve;
//std::default_random_engine generator;

typedef struct {
    std::string inputFile;
    std::string queryFile;
    std::string outputFile;
    uint16_t k;
    uint16_t M;
    uint16_t L;
    uint16_t probes;
    uint16_t N;
    double R;
    std::string algorithm;
    std::string metric;
    double delta;
}userInput;

typedef struct{
    std::string inputFile;
    std::string configFile;
    std::string outputFile;
    bool complete;
    std::string update;
    std::string assignment;
    bool silhouette;
}clusterInput;

//#endif
//#ifndef INPUT_HANDLE_HEADERS
//#define INPUT_HANDLE_HEADERS
#pragma once
#include "../headers.hpp"
#include "./LSH/hash.hpp"


int get_user_input_commant_promt(userInput *& usrIn, const int argc, char **argv);
int user_input_cluster(clusterInput *& usrIn, const int argc, char **argv);
std::vector<Point *> * parse_file(std::string fileName);
std::vector<uint32_t> * getConfigFile(std::string fileName);
void output(clusterInput* usrIn, std::unordered_map<Point*, std::unordered_set<Point*>*>* clusters, std::chrono::duration<double> duration, std::vector<double> s);
std::vector<Curve *> * add_Xaxis(std::vector<Curve *> *CurvesVector);

//#endif
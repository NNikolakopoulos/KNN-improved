#pragma once
#include "../../headers.hpp"
#include "reverse.hpp"

void lloyd(std::vector<Point*>*, std::unordered_set<Point*>*, std::unordered_map<Point*, std::unordered_set<Point*>*>*);
void lloyd_curves(std::vector<Point*>*, std::unordered_set<Point*>*, std::unordered_map<Point*, std::unordered_set<Point*>*>*);
void assignment(std::vector<Point*>*, std::unordered_set<Point*>*, std::unordered_map<Point*, std::unordered_set<Point*>*>*, std::string, userInput*, reverseAssignment*);
void assignment_curves(std::vector<Point*>*, std::unordered_set<Point*>*, std::unordered_map<Point*, std::unordered_set<Point*>*>*, std::string, userInput*, reverseAssignment*);
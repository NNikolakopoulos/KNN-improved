#include "assignment.hpp"
#include "reverse.hpp"
#include "../LSH/hash.hpp"
#include "../operations.hpp"

void lloyd(std::vector<Point*>* points, std::unordered_set<Point*>* centroids, std::unordered_map<Point*, std::unordered_set<Point*>*>* clusters)
{
    for (auto itr = centroids->begin(); itr != centroids->end(); itr++)
    {
        std::unordered_set<Point*>* cluster_points = new std::unordered_set<Point*>;
        clusters->insert({*itr, cluster_points});
    }
    for (auto itr_p = points->begin(); itr_p != points->end(); itr_p++)
    {
        auto itr_c = centroids->begin();
        auto temp = itr_c;
        double min_dist = euclideanDistance(*((*itr_p)->getVector()), *((*itr_c)->getVector()));    //distance of point from first centroid
        double cur_dist;
        itr_c++;
        while (itr_c != centroids->end())               //for every centroid
        {
            cur_dist = euclideanDistance(*((*itr_p)->getVector()), *((*itr_c)->getVector()));
            if (cur_dist < min_dist)                    //new minimum distance
            {
                min_dist = cur_dist;
                temp = itr_c;
            }
            itr_c++;
        }
        auto itr = clusters->find(*temp);
        itr->second->insert(*itr_p);
    }
    centroids->clear();
}


void lloyd_curves(std::vector<Point*>* points, std::unordered_set<Point*>* centroids, std::unordered_map<Point*, std::unordered_set<Point*>*>* clusters)
{
    for (auto itr = centroids->begin(); itr != centroids->end(); itr++)
    {
        std::unordered_set<Point*>* cluster_points = new std::unordered_set<Point*>;
        clusters->insert({*itr, cluster_points});
    }
    for (auto itr_p = points->begin(); itr_p != points->end(); itr_p++)
    {
        auto itr_c = centroids->begin();
        auto temp = itr_c;
        double min_dist = discrete_frechet(*itr_p, *itr_c);    //distance of point from first centroid
        double cur_dist;
        itr_c++;
        while (itr_c != centroids->end())               //for every centroid
        {
            cur_dist = discrete_frechet(*itr_p, *itr_c);
            if (cur_dist < min_dist)                    //new minimum distance
            {
                min_dist = cur_dist;
                temp = itr_c;
            }
            itr_c++;
        }
        auto itr = clusters->find(*temp);
        itr->second->insert(*itr_p);
    }
    centroids->clear();
}


void assignment(std::vector<Point*>* points, std::unordered_set<Point*>* centroids, std::unordered_map<Point*, std::unordered_set<Point*>*>* clusters, std::string method, userInput* ui, reverseAssignment* revAs)
{
    if (method == "Classic" || method == "classic")
    {
        lloyd(points, centroids, clusters);
    }
    else
    {
        clusters = revAs->assignment(clusters, centroids);
    }
}


void assignment_curves(std::vector<Point*>* points, std::unordered_set<Point*>* centroids, std::unordered_map<Point*, std::unordered_set<Point*>*>* clusters, std::string method, userInput* ui, reverseAssignment* revAs)
{
    if (method == "Classic" || method == "classic")
    {
        lloyd_curves(points, centroids, clusters);
    }
    else
    {
        clusters = revAs->assignment(clusters, centroids);
    }
}
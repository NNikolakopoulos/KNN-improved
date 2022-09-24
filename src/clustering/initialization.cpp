#include "../LSH/hash.hpp"
#include "../operations.hpp"
#include "initialization.hpp"

using namespace std;

unordered_set<Point*>* initialization(vector<Point*> *points, const unsigned int num_of_clusters)
{
    random_device                  rand_dev;
    mt19937                        generator(rand_dev());
    uniform_int_distribution<int>  distr(1, points->size());

    int cur_point = distr(generator);

    unordered_set<Point*> *centroids = new unordered_set<Point*>;
    unordered_set<Point*>::iterator itr_c;
    centroids->insert((*points)[cur_point-1]);
    
    for (int i = 1; i < num_of_clusters; i++)               //until num_of_clusters centroids are chosen
    {
        vector<double> distances;
        for (auto itr_p = points->begin(); itr_p != points->end(); itr_p++)         //for every point
        {
            if ( centroids->find(*itr_p) == centroids->end() )                      //if point is not a centroid
            {
                itr_c = centroids->begin();
                double min_dist = euclideanDistance(*((*itr_p)->getVector()), *((*itr_c)->getVector()));    //distance of point from first centroid
                double cur_dist;
                distances.push_back(min_dist);
                itr_c++;
                while (itr_c != centroids->end())               //for every centroid
                {
                    cur_dist = euclideanDistance(*((*itr_p)->getVector()), *((*itr_c)->getVector()));
                    if (cur_dist < min_dist)                    //new minimum distance
                    {
                        distances.pop_back();
                        distances.push_back(cur_dist);
                        min_dist = cur_dist;
                    }
                    itr_c++;
                }
            }
            else
            {
                distances.push_back(0);
            }
        }
        vector<double> part_sums(distances.size());
        part_sums[0] = distances[0] * distances[0];             //compute partial sums
        for (int j = 1; j < points->size(); j++)
        {
            part_sums[j] = part_sums[j-1] + distances[j] * distances[j];                    //we might need to normalize
        }
        uniform_real_distribution<double> distr(0, part_sums[part_sums.size()-1]);
        double num = distr(generator);
        int new_center = search_num_range(part_sums, num);           //get index to make a point new centroid
        centroids->insert((*points)[new_center]);
    }
    return centroids;
}


unordered_set<Point*>* initialization_curves(vector<Point*> *points, const unsigned int num_of_clusters)
{
    random_device                  rand_dev;
    mt19937                        generator(rand_dev());
    uniform_int_distribution<int>  distr(1, points->size());

    int cur_point = distr(generator);

    unordered_set<Point*> *centroids = new unordered_set<Point*>;
    unordered_set<Point*>::iterator itr_c;
    centroids->insert((*points)[cur_point-1]);
    
    for (int i = 1; i < num_of_clusters; i++)               //until num_of_clusters centroids are chosen
    {
        vector<double> distances;
        for (auto itr_p = points->begin(); itr_p != points->end(); itr_p++)         //for every point
        {
            if ( centroids->find(*itr_p) == centroids->end() )                      //if point is not a centroid
            {
                itr_c = centroids->begin();
                double min_dist = discrete_frechet(*itr_p, *itr_c);    //distance of point from first centroid
                double cur_dist;
                distances.push_back(min_dist);
                itr_c++;
                while (itr_c != centroids->end())               //for every centroid
                {
                    cur_dist = discrete_frechet(*itr_p, *itr_c);
                    if (cur_dist < min_dist)                    //new minimum distance
                    {
                        distances.pop_back();
                        distances.push_back(cur_dist);
                        min_dist = cur_dist;
                    }
                    itr_c++;
                }
            }
            else
            {
                distances.push_back(0);
            }
        }
        vector<double> part_sums(distances.size());
        part_sums[0] = distances[0] * distances[0];             //compute partial sums
        for (int j = 1; j < points->size(); j++)
        {
            part_sums[j] = part_sums[j-1] + distances[j] * distances[j];                    //we might need to normalize
        }
        uniform_real_distribution<double> distr(0, part_sums[part_sums.size()-1]);
        double num = distr(generator);
        int new_center = search_num_range(part_sums, num);           //get index to make a point new centroid
        centroids->insert((*points)[new_center]);
    }
    return centroids;
}


bool clusters_changed(std::list<Curve> & old_centroids, std::unordered_map<Point *, std::unordered_set<Point *> *> & clusters, const int type);


bool clusters_changed(std::list<Curve> & old_centroids, std::unordered_map<Point *, std::unordered_set<Point *> *> & clusters, const int type, const double StartRadius) {
    auto itr1 = clusters.begin();
    double total_dist = 0;
    for(auto itr2 = old_centroids.begin(); itr2 != old_centroids.end(); itr2++) {
        if( (type == TYPE_VECTOR) || (type == TYPE_CUBE) )
            total_dist += euclideanDistance( *((*itr1).first->getVector()), *((*itr2).getVector()));
        else if(type == TYPE_DISCRETE)
            total_dist = discrete_frechet((*itr1).first, &(*itr2));
        
        itr1++;
    }

    double totalAverage = total_dist / old_centroids.size(); // get average difference of distance between old and new centroids
    if( totalAverage < (StartRadius/2) )
        return false;
    else
        return true;
}

void deleteClusters(std::unordered_map<Point*, std::unordered_set<Point*>*>* clusters) {

    for(auto&& item: *clusters) {
        delete item.second;
        if (item.first->getItemId() == to_string(0))
            delete item.first;
    }
    delete clusters;
}


std::vector<double> silhouette(unordered_map<Point*, unordered_set<Point*>*>* clusters, const int typeFrechet) {

    uint32_t counter = 0;
    double s_avg = 0;
    std::vector<double> s;
    for (auto itr = clusters->begin(); itr != clusters->end(); itr++) {
        double s_cluster = 0;
        if(itr->second != nullptr) {
            for( auto it1 = itr->second->begin(); it1 != itr->second->end(); ++it1) {

                // compute a_i
                double a_i = 0.0;
                for( auto it2 = itr->second->begin(); it2 != itr->second->end(); ++it2) {

                    if(it1 == it2 )
                        continue;
                    
                    if(typeFrechet == TYPE_DISCRETE)
                        a_i += discrete_frechet( (*it1), (*it2));
                    else
                        a_i += euclideanDistance( *((*it1)->getVector()),  *((*it2)->getVector()));
                }
                a_i /= itr->second->size();

                // compute b_i


                std::unordered_set<Point*> * closest_clust_it; 
                double minDist = DBL_MAX;
                // compute closest cluster to i-point
                for(auto centroid_it = clusters->begin();  centroid_it != clusters->end(); ++centroid_it) {

                    if(centroid_it == itr )
                        continue;
                    double currDist;
                    if(typeFrechet == TYPE_DISCRETE)
                        currDist = discrete_frechet( centroid_it->first, itr->first);
                    else
                        currDist = euclideanDistance( *(centroid_it->first->getVector()), *(itr->first->getVector()) ); 
                    if(currDist < minDist) {
                        minDist = currDist;
                        closest_clust_it = centroid_it->second;
                    }
                }

                double b_i = 0.0;
                for( auto it2 = closest_clust_it->begin(); it2 != closest_clust_it->end(); ++it2) {

                    if(it1 == it2 )
                        continue;
                    
                    if(typeFrechet == TYPE_DISCRETE)
                        b_i += discrete_frechet( (*it1), (*it2));
                    else
                        b_i += euclideanDistance( *((*it1)->getVector()),  *((*it2)->getVector()));
                }
                b_i /= itr->second->size();

                //compute s_i

                double s_i = ( b_i - a_i ) / (std::max(a_i,b_i));


                s_cluster += s_i;
                //s_avg += s_i;
                //counter++;
            }
            s_cluster /= itr->second->size();
        }
        else 
            s_cluster = 0;
        s_avg += s_cluster;
        counter++;
        s.push_back(s_cluster);
        
    }
    s_avg /= static_cast<double>(counter);
    s.push_back(s_avg);
    return s;
}

void deletePointsVector(std::vector<Point *> * vector, std::string method) {
    if (method == "Classic" || method == "classic")
    {
        for(auto&& item: *vector) {
            delete item;
        }
    }
    delete vector;
}
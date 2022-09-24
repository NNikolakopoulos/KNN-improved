
#include "reverse.hpp"


reverseAssignment::reverseAssignment(std::vector<Curve *> *CurvesVector, userInput *ui, std::unordered_set<Curve *> *centroids, const int method) 
    : lsh((method==TYPE_VECTOR) ? (new LSH(CurvesVector, ui, TYPE_VECTOR)) : nullptr ) ,
      lsh_frechet((method==TYPE_DISCRETE) ? (new LSH_discrete_frechet(CurvesVector, ui, compute_delta())) : nullptr),
      cube((method==TYPE_CUBE) ? (new RandomizedProjection(CurvesVector, ui)) : nullptr ) ,  
      centroids(centroids) , 
      method(method)
    {}

reverseAssignment::~reverseAssignment() {
    if(method == TYPE_VECTOR) {
        delete lsh;
    }
    else if(method == TYPE_CUBE)
        delete cube;
    else if(method == TYPE_DISCRETE)
        delete lsh_frechet;
}


double reverseAssignment::calculateStartRadius() const {
 
    double minDist = DBL_MAX;

    for (auto itr1 = centroids->begin(); itr1 != centroids->end(); itr1++) {
        for (auto itr2 = itr1; itr2 != centroids->end(); itr2++) {
            if(itr2 == itr1)
                continue;
            double currDist;
            if( (this->method == TYPE_VECTOR) || (this->method == TYPE_CUBE))
                currDist = euclideanDistance( *(*itr1)->getVector(), *(*itr2)->getVector());
            else if(this->method == TYPE_DISCRETE)
                currDist = discrete_frechet((*itr1), (*itr2));
            if( currDist < minDist )
                minDist = currDist;
        }
    }

    return minDist/2;
}


std::unordered_map<Point *, std::unordered_set<Point *> *> *reverseAssignment::assignment(std::unordered_map<Point *, std::unordered_set<Point *> *> *clusters, 
 std::unordered_set<Point *> *centroids) {

    double R = this->calculateStartRadius();
    double StartRadius = R;
    bool addBarrier = true;
    //std::list<Curve > old_centroids;
    int clusters_changed = clusters->size() + 1;

    do {                                  // while at least half of the clusters have significant changes
        clusters_changed = 0;
        for( auto it = centroids->begin(); it != centroids->end(); ++it ) {             // for each centroid

            //old_centroids.push_back(**it);
            Point *centroid = *it;
            auto it_cluster = (*clusters).find(centroid);                             // get the cluster of this centroid
            if(it_cluster == (*clusters).end())                                        // if we didnt find cluster for this centroid, then initialize
                clusters->insert({centroid , new std::unordered_set<Point *>});

            std::unordered_set<Point *> * cluster = (*(*clusters).find(centroid)).second;
            if(this->method == TYPE_VECTOR)
                clusters_changed = (this->lsh->rangedSearchCentroid(centroid, cluster, R, addBarrier)) ? (clusters_changed + 1) : clusters_changed;
            else if(this->method == TYPE_DISCRETE)
                clusters_changed = (this->lsh_frechet->rangedSearchCentroidFrechet(centroid, cluster, R, addBarrier)) ? (clusters_changed + 1) : clusters_changed;
            else if(this->method == TYPE_CUBE)
                clusters_changed = (this->cube->rangedSearchCentroid(centroid, cluster, R, addBarrier)) ? (clusters_changed + 1) : clusters_changed; 
            /*if(this->method == TYPE_VECTOR)
                this->lsh->rangedSearchCentroid(centroid, cluster, R, addBarrier);
            else if(this->method == TYPE_DISCRETE)
                this->lsh_frechet->rangedSearchCentroidFrechet(centroid, cluster, R, addBarrier);
            else if(this->method == TYPE_CUBE)
                this->cube->rangedSearchCentroid(centroid, cluster, R, addBarrier); */
            addBarrier = false;
        }
        R *= 2;                                                                         // increase radius
    }while(clusters_changed > (centroids->size()/2));

    if(this->method == TYPE_VECTOR)
        this->lsh->unassignedPoints(clusters, centroids);
    else if(this->method == TYPE_CUBE)
        this->cube->unassignedPoints(clusters, centroids);
    else if(this->method == TYPE_DISCRETE)
        this->lsh_frechet->unassignedPoints(clusters, centroids);

    return clusters;
}



#include "hyperc.hpp"




HyperCube::HyperCube(const uint32_t tableSize) 
    : array(new std::vector<std::list<Point *>*>(tableSize, NULL)), size(tableSize) {}

HyperCube::~HyperCube() {
     for(int i=0; i<size; i++) {
        std::list<Point *> *list = (*array)[i];

        if(list == NULL)
            continue;

        for(auto&& item : *list) {
           delete item;
        }

        delete list;
    }
    delete array;
}



const int HyperCube::insert(Point *pointPtr, const uint64_t key) {
    if((*array)[key] != NULL) {
        std::list<Point *>::iterator it;
        for( it = (*array)[key]->begin(); it != (*array)[key]->end(); ++it ) {
            if((*it)->getItemId() == pointPtr->getItemId())                // dont insert duplicates
                return 0;
        }
    }
    else
        (*array)[key] = new std::list<Point *>;

    (*array)[key]->push_back(pointPtr);                                    // now insert the point 
    ids.insert({pointPtr->getItemId(),key});
    return 1;
}   

std::multimap<double, std::pair<Point *, std::chrono::duration<double>> > & HyperCube::bruteForce(Point *pointPtr, const uint32_t K, std::multimap<double, std::pair<Point *, std::chrono::duration<double>> > & map) const {

    double KminDist, currDist;
    Point *minPoint;
    std::list<Point *>::iterator it;
    std::list<Point *> *list;

    KminDist = DBL_MAX;

    for(int i=0; i<size; ++i) {
        list = (*array)[i];
        if(list != NULL) {
            auto start = std::chrono::high_resolution_clock::now();
            for(it = list->begin(); it != list->end(); ++it) {               // for each point in bucket
                currDist = euclideanDistance( *((*it)->getVector()) , *(pointPtr->getVector()) );
                if(currDist < KminDist) {
                    auto stop = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> duration = stop - start;
                    if(K <= map.size())  {                           // if we have K neighbors, then delete last(most distant) {
                        map.erase(std::prev(map.end()));            
                    }
                    map.insert( { currDist, std::pair<Point *, std::chrono::duration<double>> (*it, duration) });
                    KminDist = map.rbegin()->first;
                    minPoint = map.rbegin()->second.first;

                }
            }
        }
    }
    return map;
}


void HyperCube::KNN(Point *pointPtr, const u_int32_t key, const uint32_t K,std::multimap<double,std::pair<Point *, std::chrono::duration<double>>> & map, userInput *ui) const {

    uint32_t newKey, countM, countProbes;
    double KminDist, currDist;
    Point *minPoint;
    std::list<Point *>::iterator it;
    std::list<Point *> *list;
    std::queue<uint32_t> queue;
    std::unordered_set<uint32_t> expl_set_hamDist;

    queue.push(key);
    queue = hamDist1(queue, key, expl_set_hamDist, ui->k);

    countM = 0;
    countProbes = 0;
    KminDist = DBL_MAX;

    while(countProbes < ui->probes && countM < ui->M) {

        if(queue.empty()) {
            break;

        }

        newKey = queue.front();
        queue.pop();
        
        list = (*array)[newKey];

        if(list != NULL ) {
            for(it = list->begin(); it != list->end(); ++it) {                  // for each point in bucket
            auto start = std::chrono::high_resolution_clock::now();
                if(countM >= ui->M)
                    break;
                if(ids.find((*it)->getItemId())->second == key) {
                    currDist = euclideanDistance( *((*it)->getVector()) , *(pointPtr->getVector()) );
                    if(currDist < KminDist) {
                        auto stop = std::chrono::high_resolution_clock::now();
                        std::chrono::duration<double> duration = stop - start;
                        if(K <= map.size())                                          // if we have K neighbors, then delete last(most distant)
                            map.erase(std::prev(map.end()));            
                        map.insert({currDist , std::pair<Point *, std::chrono::duration<double>> (*it, duration) } );
                        KminDist = map.rbegin()->first;
                        minPoint = map.rbegin()->second.first;

                    }
                }
                countM++;
            }
        }
        countProbes++;
        queue = hamDist1(queue, newKey, expl_set_hamDist, ui->k);
    }
}


// range search but inserts barrier and stops searching at barrier of every buvket. Returns number of points inserted in cluster
int HyperCube::rangedSearchAmplified(Point *pointPtr, const uint64_t key, std::unordered_set<Point *> *point_set, userInput *ui, const double R, const bool addBarrier) 
{
    int changed_cluster_flag;
    uint32_t newKey, countM, countProbes;
    double KminDist, currDist;
    Point *minPoint;
    std::list<Point *>::iterator it, iter;
    std::list<Point *> *list;
    std::queue<uint32_t> queue;
    std::unordered_set<uint32_t> expl_set_hamDist;

    queue.push(key);
    queue = hamDist1(queue, key, expl_set_hamDist, ui->k);

    countM = 0;
    countProbes = 0;
    changed_cluster_flag = 0;
    KminDist = DBL_MAX;

    while(countProbes < ui->probes && countM < ui->M) {

        if(queue.empty()) {
            break;
        }

        newKey = queue.front();
        queue.pop();
        
        list = (*array)[newKey];

        if(list != NULL ) {
            if(addBarrier)                                                              // inside this bucket, we will move points with dist < R after this barrier  
                list->push_back(new Point("-1",nullptr));

            for(it = list->begin(); it != list->end(); ++it) {                  // for each point in bucket

                if(countM >= ui->M)
                    break;
                    
                if( (*it)->getVector() == nullptr )                                     // if we found barrier then stop.
                    return changed_cluster_flag;
                if(point_set->find((*it)) == point_set->end()) {

                    currDist = euclideanDistance( *((*it)->getVector()) , *(pointPtr->getVector()) );

                    if(currDist < R) {
                        changed_cluster_flag++;
                        point_set->insert(*it);                                     // insert this point to the given centroid cluster
                        list->splice(list->end() , *list, it);                    // move this item to the end of the bucket ( after the barrier )
                    }
                }

                countM++;
            }
        }
        countProbes++;
        queue = hamDist1(queue, newKey, expl_set_hamDist, ui->k);
    }
    return changed_cluster_flag;
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Randomized Projection ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


RandomizedProjection::RandomizedProjection(std::vector<Point *> * PointsVector, userInput *ui) 
    : ui(ui) , hyperCube(new HyperCube(pow2(ui->k)))
{
    uint32_t tableSize = PointsVector->size();
    Point *point = (*PointsVector)[0];
    uint32_t dim = point->getVector()->size();

    //uint16_t w = 747;                           // window value
    uint16_t w = 200; 
    uint32_t HTsize = pow2(ui->k);
    this->hashFun = new HashFs(ui->k, w , HTsize, dim);
    int duplicates=0;                           // now insert all points to hash tables

    for(int j=0; j < tableSize; j++) {      // for all points
        point = (*PointsVector)[j];
        uint64_t key = this->hashFun->amplified_hash(point);
        if( this->hyperCube->insert(point,key) == 0)  
            duplicates++;
    }
    //delete PointsVector;
}

RandomizedProjection::~RandomizedProjection() {
    delete ui;
    delete hashFun;
    delete hyperCube;
}


std::multimap<double, std::pair<Point *, std::chrono::duration<double>>> RandomizedProjection::queryKNN(Point *pointPtr, const uint32_t K) const {
    std::string itemId;
    std::tuple<Point *, double> tuple;
    std::multimap<double, std::pair<Point *, std::chrono::duration<double>>> map;
    double minDist = DBL_MAX, currDist;

    uint64_t key = this->hashFun->amplified_hash(pointPtr);      // get hashed key from amplified function for query point

    hyperCube->KNN(pointPtr, key, K, map, ui);              // get the nearest K-points and its distances via map

    return map;
}



bool RandomizedProjection::queryFile(std::string queryFile) const {

    std::ofstream outputFile(ui->outputFile);
    std::vector<Point *> *queryPoints;
    queryPoints = parse_file(queryFile);
    if(queryPoints == nullptr)
        return false;
    double MAF = 0.0001;
    for(int i=0; i<queryPoints->size(); i++) {
        Point *point =(*queryPoints)[i];
        //std::tuple<std::string, double>  NN = queryNN(point);
        outputFile << "\nQuery: "+point->getItemId() << '\n'; //<< "Nearest Neighbor-1: "+ std::get<0>(NN) << "\nDistance LSH: "+  std::to_string(std::get<1>(NN))<<"\n\n";
        std::multimap<double,std::pair<Point *, std::chrono::duration<double>>> map = queryKNN(point, ui->N);
        std::multimap<double, std::pair<Point *, std::chrono::duration<double>>> bruteForceMap;
        bruteForceMap = hyperCube->bruteForce(point, ui->N, bruteForceMap);
        std::multimap<double, std::pair<Point *, std::chrono::duration<double>>>::iterator it;
        std::multimap<double, std::pair<Point *, std::chrono::duration<double>>>::iterator bruteForceIterator = bruteForceMap.begin();
        int count=0;
        for(it = map.begin(); it != map.end(); ++it) {
            MAF = (it->first/bruteForceIterator->first  > MAF) ? it->first/bruteForceIterator->first : MAF;
            outputFile << "Nearest neighbor-"<< ++count << ": " << it->second.first->getItemId()<< "\nDistanceLSH: "<<  it->first << "\nDistanceTrue: " << bruteForceIterator->first <<'\n' ;
            outputFile << "tLSH: "<<  it->second.second.count() << "\ntTrue: " << bruteForceIterator->second.second.count() <<'\n' ;
            bruteForceIterator++;
        }
        delete point;
    }
    delete queryPoints;
    outputFile << "MAF: "<< MAF << '\n';
    outputFile.close();
    return true;
}

int RandomizedProjection::rangedSearchCentroid(Point *centroid, std::unordered_set<Point *> *cluster, const double R, const bool addBarrier) {

    std::unordered_set<std::string> explored_set;
    double minDist = DBL_MAX;
    int points_inserted_into_clusters = 0;


    uint64_t key = this->hashFun->amplified_hash(centroid);      // get hashed key from amplified function for query point
    points_inserted_into_clusters += hyperCube->rangedSearchAmplified(centroid, key, cluster, ui,  R, addBarrier);              // get the nearest K-points and its distances via map

    return points_inserted_into_clusters;
    /*if(points_inserted_into_clusters < ( cluster->size()/10 ) )         // if only little number of new points are inserted into cluster
        return false;                                                   // return false, to stop the assignment
    else
        return true;  */
}

// assign all point that they aren't in any cluster with bruteforce, delete barriers and if any point is in >=2 clusters, then assign to the closer one
void RandomizedProjection::unassignedPoints(std::unordered_map<Point *, std::unordered_set<Point *> *> *clusters, std::unordered_set<Point *> *centroids) {

    std::vector<std::list<Point *>*> *points =hyperCube->getArray();
    std::list<Point *> concurrentCentroids;

    for(int i=0; i<points->size(); i++) {                                                                                   // for every bucket 

        std::list<Point *> *list = points->at(i);   
        if(list == nullptr)
            continue;
        auto list_it = list->begin();
        while( list_it != list->end() ) {                                              // for every point inside this bucket  

            if((*list_it)->getVector() == nullptr)  {                                   // if its barrier then delete it
                list->erase(list_it++);                       
                continue;
            }
            concurrentCentroids.clear();
            for(auto centroid_it = centroids->begin(); centroid_it != centroids->end(); ++centroid_it) {                    // for every centroid
                Point *centroid = *centroid_it;
                auto cluster_it =  clusters->find(centroid);
                 
                // if this cluster exists and the point exists in this cluster
                if(cluster_it != clusters->end() && cluster_it->second->find(*list_it) != cluster_it->second->end() )
                    concurrentCentroids.push_back(centroid);
            }
            if(concurrentCentroids.size() == 0 )
                assignPoint(clusters, *list_it, TYPE_VECTOR);
            else if(concurrentCentroids.size() >= 2) 
                assignPointCloserCluster(clusters, concurrentCentroids, *list_it, TYPE_VECTOR);

            list_it++;
        }
    }
}

#include "hash.hpp"


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   POINT   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


Point::Point(std::tuple<std::string, std::vector<double> *> tuple) 
    :vect(std::get<1>(tuple)), itemId(std::get<0>(tuple)) {}

Point::Point(std::string id, std::vector<double> *vector) 
    :vect(vector), itemId(id) {}

Point::~Point() {
    delete vect;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   HASHTABLE   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */



HashTable::HashTable(const uint32_t tableSize) 
    : array(new std::vector<std::list<Point *>*>(tableSize, NULL)), size(tableSize) {}

HashTable::~HashTable() {

    for(int i=0; i<size; i++) {
        std::list<Point *> *list = (*array)[i];

        if(list == NULL)
            continue;

        /*
        for(auto&& item : *list) {
           delete item;
        }
        for( auto it = list->begin(); it != list->end(); ++it) {
            (*it)->~Point();
            //delete *it;
        }*/

        delete list;
    }
    delete array;
}

void HashTable::deletePoints() {

    for(int i=0; i<size; i++) {
        std::list<Point *> *list = (*array)[i];
        
        if(list == NULL)
            continue;
        for(auto&& item : *list) {
            delete item;
        }
    }
}

const int HashTable::insert(Point *pointPtr, const uint64_t key) {
    uint64_t safeKey = mod(key,size);
    if((*array)[safeKey] != NULL) {
        std::list<Point *>::iterator it;
        //for( it = (*array)[safeKey]->begin(); it != (*array)[safeKey]->end(); ++it ) {
        //    if((*it)->getItemId() == pointPtr->getItemId())                // dont insert duplicates
        //        return 0;
        //}
    }
    else
        (*array)[safeKey] = new std::list<Point *>;

    (*array)[safeKey]->push_back(pointPtr);                                    // now insert the point 
    ids.insert({pointPtr->getItemId(),key});
    return 1;
}   



std::multimap<double, std::pair<Point *, std::chrono::duration<double>> > & HashTable::bruteForce(Point *pointPtr, const uint32_t K,std::multimap<double, std::pair<Point *, std::chrono::duration<double>> > & map, int typeFrechet) const {

    double KminDist, currDist;
    Point *minPoint;
    std::list<Point *>::iterator it;
    std::list<Point *> *list;

    KminDist = DBL_MAX;
    
    for(int i=0; i<size; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        list = (*array)[i];
        if(list != NULL) {
            for(it = list->begin(); it != list->end(); ++it) {               // for each point in bucket
                if(typeFrechet == 0)
                    currDist = euclideanDistance( *((*it)->getVector()) , *(pointPtr->getVector()) );
                else if(typeFrechet == TYPE_DISCRETE)
                    currDist = discrete_frechet( (*it) , pointPtr );
                else if(typeFrechet == TYPE_CONTINUOUS) {
                    fred::Curve * fred_curve1 = convertCurve((*it)->getVector());
                    fred::Curve * fred_curve2 = convertCurve(pointPtr->getVector());
                    currDist = fred::Frechet::Continuous::distance(*fred_curve1, *fred_curve2).value;
                    delete fred_curve1;
                    delete fred_curve2;
                }
                if(currDist < KminDist) {
                    auto stop = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> duration = stop - start;
                    if(K <= map.size())  {                           // if we have K neighbors, then delete last(most distant) {
                        map.erase(std::prev(map.end()));            
                    }
                    map.insert( { currDist, std::pair<Point *, std::chrono::duration<double>> (*it, duration) });
                    KminDist = map.rbegin()->first;
                    minPoint = map.rbegin()->second.first;
                    auto start = std::chrono::high_resolution_clock::now();
                }
                auto start = std::chrono::high_resolution_clock::now();
            }
        }
    }
    return map;
}


void HashTable::KNN(Point *pointPtr, const u_int64_t key, const uint32_t K,std::multimap<double, std::pair<Point *, std::chrono::duration<double>> > & map, std::unordered_set<std::string > & explored_set, int typeFrechet) const {

    double KminDist, currDist;
    Point *minPoint;
    std::list<Point *>::iterator it;
    std::list<Point *> *list = (*array)[mod(key,size)];
    if(list == NULL)            // if this bucket has no points then return
        return;
    //KminDist = (!map.size()) ? DBL_MAX : map.rbegin()->first;
    if(map.size() == 0)         // if we havent found any
        KminDist = DBL_MAX;
    else
        KminDist = map.rbegin()->first; 
    for(it = list->begin(); it != list->end(); ++it) {               // for each point in bucket
        auto start = std::chrono::high_resolution_clock::now();
        if( (ids.find((*it)->getItemId())->second == key) || (typeFrechet == TYPE_CONTINUOUS) ) {       // check points with the same ID(p) = (r1*h1(p)+...+ri*hi(p)) mod M
            if(explored_set.find((*it)->getItemId()) == explored_set.end()) {       // if this point is not explored yet from other hashtables, then find eucl.dist.
                if(typeFrechet == 0)
                    currDist = euclideanDistance( *((*it)->getVector()) , *(pointPtr->getVector()) );
                else if(typeFrechet == TYPE_DISCRETE)
                    currDist = discrete_frechet( (*it) , pointPtr );
                else if(typeFrechet == TYPE_CONTINUOUS) {
                    fred::Curve * fred_curve1 = convertCurve((*it)->getVector());
                    fred::Curve * fred_curve2 = convertCurve(pointPtr->getVector());
                    currDist = fred::Frechet::Continuous::distance(*fred_curve1, *fred_curve2).value;
                    delete fred_curve1;
                    delete fred_curve2; 
                    std::cout << currDist << '\n';
                }
                if(currDist < KminDist) {
                    auto stop = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> duration = stop - start;
                    if(K <= map.size())                             // if we have K neighbors, then delete last(most distant)
                        map.erase(std::prev(map.end()));            
                    map.insert( { currDist, std::pair<Point *, std::chrono::duration<double>>(*it,duration) } );
                    KminDist = map.rbegin()->first;
                    minPoint = map.rbegin()->second.first;
                }
                explored_set.insert((*it)->getItemId());
            }
        }
    }
}


int HashTable::rangedSearchAmplified(Point *pointPtr, const uint64_t key, const double R, std::unordered_set<Point *> *point_set, 
            std::unordered_set<std::string> & explored_set, const bool addBarrier, const int typeFrechet)
{
    int changed_cluster_flag = 0;
    double currDist;
    std::list<Point *>::iterator it, iter;
    std::list<Point *> *list = (*array)[mod(key,size)];

    if(list == NULL)                                                            // if this bucket has no points then return
        return changed_cluster_flag;

    if(addBarrier)                                                              // inside this bucket, we will move points with dist < R after this barrier  
        list->push_back(new Point("-1",nullptr));
    for(it = list->begin(); it != list->end(); ++it) {                          // for each point in bucket
        if( (*it)->getVector() == nullptr )                                     // if we found barrier then stop.
            return changed_cluster_flag;

        if( (ids.find((*it)->getItemId())->second == key) || (TYPE_DISCRETE == typeFrechet)) {                      // check points with the same ID(p) = (r1*h1(p)+...+ri*hi(p)) mod M
            if(explored_set.find((*it)->getItemId()) == explored_set.end()) {       // if this point is not explored yet from other hashtables, then find eucl.dist.
                if(point_set->find((*it)) == point_set->end()) {

                    std::string itemId = (*it)->getItemId();
                    if(typeFrechet == TYPE_VECTOR)
                        currDist = euclideanDistance( *((*it)->getVector()) , *(pointPtr->getVector()) );
                    else if(typeFrechet == TYPE_DISCRETE)
                        currDist = discrete_frechet( (*it) , pointPtr );
                    if(currDist < R) {
                        changed_cluster_flag++;
                        point_set->insert(*it);                                     // insert this point to the given centroid cluster
                        list->splice(list->end() , *list, it);                    // move this item to the end of the bucket ( after the barrier )
                    }
                    explored_set.insert(itemId);
                }
            }
        }
    }
    return changed_cluster_flag;
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   LSH   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

LSH::LSH(std::vector<Point *> * PointsVector, userInput *ui, int typeFrechet) : typeFrechet(typeFrechet), ui(ui)
{
    uint32_t tableSize = PointsVector->size();
    Curve *curve =  (*PointsVector)[0];
    uint32_t dim = curve->getVector()->size();

    uint64_t M = UINT32_MAX - 5;                // M value
    uint16_t w = 100;                             // window value
    uint32_t HTsize = tableSize/8;

    this->hashFun = new HashFunctions(ui->k, w , M, ui->L,HTsize, dim);
    int duplicates=0;                           // now insert all points to hash tables

    // we have only 1 hash table per LSH
    for(int i=0; i<ui->L; i++)
        HTarray.push_back(new HashTable(HTsize));

    for(int i=0; i<ui->L; i++)
        for(int j=0; j < tableSize; j++) {      // for all points
            curve = (*PointsVector)[j];

            uint64_t key = this->hashFun->amplified_hash(curve, i);


            if( HTarray[i]->insert(curve,key) == 0 && i==0)   // count duplicates only 1 time
                duplicates++;                
        }
    //delete PointsVector;
}

// CONSTRUCTOR FOR DISCRETE FRESHET
LSH::LSH(std::vector<Point *> * PointsVector, userInput *ui, const double delta, std::pair<double,double> *t, int typeFrechet) : typeFrechet(typeFrechet), ui(ui)
{
    uint32_t tableSize = PointsVector->size();
    Curve *curve =  (*PointsVector)[0];
    uint32_t dim = curve->getVector()->size();

    uint64_t M = UINT32_MAX - 5;                // M value
    uint16_t w = 150;                             // window value
    uint32_t HTsize = tableSize/8;

    this->hashFun = new HashFunctions(ui->k, w , M, ui->L,HTsize, dim);
    int duplicates=0;                           // now insert all points to hash tables

    // we have only 1 hash table per LSH
    HTarray.push_back(new HashTable(HTsize));

    for(int j=0; j < tableSize; j++) {      // for all points
        curve = (*PointsVector)[j];

        std::vector<std::pair<double,double>*> *snapped_curve = snap(*curve, delta, *t);
        std::vector<double> * concat_curve = getConcatCurve(snapped_curve);
        uint64_t key = this->hashFun->amplified_hash(concat_curve, 0);
        delete concat_curve;

        if( HTarray[0]->insert(curve,key) == 0 )   // count duplicates only 1 time
            duplicates++;    
    }
    //delete PointsVector;
}


// CONSTRUCTOR FOR CONTINUOUS FRECHET
LSH::LSH(std::vector<Curve *> * CurvesVector, userInput *ui, const double delta, double t , double epsilon,int typeFrechet) : typeFrechet(typeFrechet), ui(ui)
{
    uint32_t tableSize = CurvesVector->size();
    Curve *curve =  (*CurvesVector)[0];
    uint32_t dim = curve->getVector()->size();

    uint64_t M = UINT32_MAX - 5;                // M value
    uint16_t w = 100;                             // window value
    uint32_t HTsize = tableSize/8;

    this->hashFun = new HashFunctions(ui->k, w , M, ui->L,HTsize, dim);
    int duplicates=0;                           // now insert all points to hash tables

    // we have only 1 hash table 
    HTarray.push_back(new HashTable(HTsize));

    for(int j=0; j < tableSize; j++) {      // for all points
        curve = (*CurvesVector)[j];

        std::vector<double> * filtered_curve = filter(*curve, epsilon);
        filtered_curve = snap_continuous(filtered_curve, delta, t);
        filtered_curve = minima_maxima(filtered_curve);

        uint64_t key;
        key = this->hashFun->amplified_hash(filtered_curve, 0);


        if( HTarray[0]->insert(curve,key) == 0)   // count duplicates only 1 time
            duplicates++;                
    }
    //delete CurvesVector;
}

LSH::~LSH() {

    //HTarray[0]->deletePoints();
    for(int i=0; i<HTarray.size(); i++) {
        delete HTarray[i];
    }

    delete hashFun;
    //delete ui;
}

std::multimap<double, std::pair<Point *, std::chrono::duration<double>> >& LSH::queryKNN(Point *pointPtr, const uint32_t K, 
  std::multimap<double, std::pair<Point *, std::chrono::duration<double>> >&map, std::unordered_set<std::string> &explored_set, const double delta, std::pair<double,double> *t, const double epsilon) const 
{

    int numHashTables = (typeFrechet != 0) ? 1 : ui->L;
    for(int i=0; i<numHashTables; i++) {
        //uint64_t key = this->hashFun->amplified_hash(pointPtr, 0);      // get hashed key from amplified function for query point
        uint64_t key;
        if(typeFrechet == TYPE_VECTOR)
            key = this->hashFun->amplified_hash(pointPtr, i);
        else if(typeFrechet == TYPE_DISCRETE) {
            std::vector<std::pair<double,double>*> *snapped_curve = snap(*pointPtr, delta, *t);
            std::vector<double> * concat_curve = getConcatCurve(snapped_curve);
            key = this->hashFun->amplified_hash(concat_curve, i);
            delete concat_curve;
        }
        else if(typeFrechet == TYPE_CONTINUOUS) {
            std::vector<double> * filtered_curve = filter(*pointPtr , epsilon);
            filtered_curve = snap_continuous(filtered_curve, delta, t->first);
            filtered_curve = minima_maxima(filtered_curve);
            key = this->hashFun->amplified_hash(filtered_curve, i);
        }
        HTarray[i]->KNN(pointPtr, key, K, map, explored_set, this->typeFrechet);              // get the nearest K-points and its distances via map
    }
    return map;
}


bool LSH::queryFile(std::string queryFile) const {

    std::ofstream outputFile(ui->outputFile);
    std::vector<Point *> *queryPoints;
    queryPoints = parse_file(queryFile);
    if(queryPoints == nullptr)      // if invalid query file is given return false
        return false;
    
    std::multimap<double, std::pair<Point *, std::chrono::duration<double>>> map;
    std::unordered_set<std::string> explored_set;
    double MAF = 0.001;
    for(int i=0; i<queryPoints->size(); i++) {
        Point *point = (*queryPoints)[i];
        //std::tuple<std::string, double>  NN = queryNN(point);
        outputFile << "Query: "+point->getItemId() << '\n'; //<< "Nearest Neighbor-1: "+ std::get<0>(NN) << "\nDistance LSH: "+  std::to_string(std::get<1>(NN))<<"\n\n";
        map = queryKNN(point, ui->N, map, explored_set, 0, nullptr, 0);
        std::multimap<double, std::pair<Point *, std::chrono::duration<double>>> bruteForceMap;
        bruteForceMap = HTarray[0]->bruteForce(point, ui->N, bruteForceMap, typeFrechet);
        std::multimap<double, std::pair<Point *, std::chrono::duration<double>>>::iterator it;
        std::multimap<double, std::pair<Point *, std::chrono::duration<double>>>::iterator bruteForceIterator = bruteForceMap.begin();
        int count=0;
        for(it = map.begin(); it != map.end(); ++it) {
            MAF = (it->first/bruteForceIterator->first  > MAF) ? it->first/bruteForceIterator->first : MAF;
            outputFile << "Nearest neighbor-"<< ++count << ": " << it->second.first->getItemId()<< "\nDistanceLSH: "<<  it->first << "\nDistanceTrue: " << bruteForceIterator->first <<'\n' ;
            outputFile << "tLSH: "<<  it->second.second.count() << "\ntTrue: " << bruteForceIterator->second.second.count() <<'\n' ;
            bruteForceIterator++;
        }
        map.clear();
        explored_set.clear();
        delete point;

    }
    delete queryPoints;
    outputFile << "MAF: "<< MAF << '\n';
    outputFile.close();
    return true;
}

int LSH::rangedSearchCentroid(Point * centroid, std::unordered_set<Point *> *cluster, const double R, const bool addBarrier) {
    
    std::unordered_set<std::string> explored_set;
    double minDist = DBL_MAX;
    int points_inserted_into_clusters = 0;

    for(int i=0; i<ui->L; i++) {

        uint64_t key = this->hashFun->amplified_hash(centroid, i);      // get hashed key from amplified function for query point
        points_inserted_into_clusters += HTarray[i]->rangedSearchAmplified(centroid, key, R, cluster, explored_set, addBarrier, TYPE_VECTOR);              // get the nearest K-points and its distances via map

    }
    return points_inserted_into_clusters;
}

int LSH::rangedSearchCentroid(Point * centroid, std::unordered_set<Point *> *cluster, const double R, const bool addBarrier,std::unordered_set<std::string> &explored_set, const double delta, std::pair<double,double>& t ) {
    
    int points_inserted_into_clusters;

    std::vector<std::pair<double,double>*> *snapped_curve = snap(*centroid, delta, t);
    std::vector<double> * concat_curve = getConcatCurve(snapped_curve);
    uint64_t key = this->hashFun->amplified_hash(concat_curve, 0);
    delete concat_curve;
    points_inserted_into_clusters = HTarray[0]->rangedSearchAmplified(centroid, key, R, cluster, explored_set, addBarrier, TYPE_DISCRETE);              // get the nearest K-points and its distances via map

    return points_inserted_into_clusters;
}


// assign all point that they aren't in any cluster with bruteforce, delete barriers and if any point is in >=2 clusters, then assign to the closer one
void LSH::unassignedPoints(std::unordered_map<Point *, std::unordered_set<Point *> *> *clusters, std::unordered_set<Point *> *centroids) {

    std::vector<std::list<Point *>*> *points = HTarray[0]->getArray();
    std::list<Point *> concurrentCentroids;

    for(int i=0; i<points->size(); i++) {                                                                                   // for every bucket 

        std::list<Point *> *list = points->at(i);   
        if(list == nullptr)
            continue;
        auto list_it = list->begin();
        while( list_it != list->end() ) {                                              // for every point inside this bucket  

            if((*list_it)->getVector() == nullptr) {                                    // if its barrier then delete it
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
                assignPoint(clusters, *list_it, this->typeFrechet);
            else if(concurrentCentroids.size() >= 2) 
                assignPointCloserCluster(clusters, concurrentCentroids, *list_it, this->typeFrechet);

            list_it++;
        }
    }
}

// insert point to the closest centroid using bruteforce
void assignPoint(std::unordered_map<Point *, std::unordered_set<Point *> *> *clusters, Point *point, const int typeFrechet) {

    double currDist, minDist;
    Point *closestCentroid;

    minDist = DBL_MAX;
    for(auto it = clusters->begin(); it != clusters->end(); ++it ) {
        Point *centroid = it->first;
        if(typeFrechet == TYPE_VECTOR)
            currDist = euclideanDistance( *(centroid->getVector()) , *(point->getVector()) );
        else if(typeFrechet == TYPE_DISCRETE)
            currDist = discrete_frechet(centroid, point);
        if(currDist < minDist) {
            minDist = currDist;
            closestCentroid = centroid;
        }
    }
    clusters->find(closestCentroid)->second->insert(point);
}

// if a point exists in more than 1 clusters, find the closest and delete it from the rest
void assignPointCloserCluster(std::unordered_map<Point *, std::unordered_set<Point *> *> *clusters, std::list<Point *> &concurrentCentroids, Point *point, const int typeFrechet) {

    double currDist, minDist;
    Point *closestCentroid;

    //calculate the closest centroid
    minDist = DBL_MAX;
    for(auto it = concurrentCentroids.begin(); it != concurrentCentroids.end(); ++it ) {
        Point *centroid = *it;
        if(typeFrechet == TYPE_VECTOR)
            currDist = euclideanDistance( *(centroid->getVector()) , *(point->getVector()) );
        else if(typeFrechet == TYPE_DISCRETE)
            currDist = discrete_frechet(centroid, point);

        if(currDist < minDist) {
            minDist = currDist;
            closestCentroid = centroid;   
        }
    }

    // delete the centroid from all clusters except one(the closest)
    auto it = concurrentCentroids.begin();
    while( it != concurrentCentroids.end()) {
        Point *centroid = *it;

        if(*it != closestCentroid)
            clusters->find(*it)->second->erase(point);
            //concurrentCentroids.erase(it++);
        it++;
    }
}


std::vector<std::pair<double,double>*> * LSH::snap(Curve &curve, const double delta, std::pair<double,double> &t) const{

    int numPoints = (curve.getVector()->size())/2;
    std::vector<std::pair<double,double>*> *snapped_curve = new std::vector<std::pair<double,double>*>(numPoints);
    std::pair<double,double>* pair;

    int duplicates = 0;

    int i=0, totalPoints = 0;
    while(totalPoints<numPoints) {
        totalPoints++;
        double x = curve.getVector()->at(i+numPoints);
        double y = curve.getVector()->at(i);
        double new_x = ((x-t.first)/delta + 0.5)*delta + t.first;
        double new_y = ((y-t.second)/delta + 0.5)*delta + t.second;

        if(i>0 && snapped_curve->at(i-1)->first == new_x && snapped_curve->at(i-1)->second == new_y) {   // check if consecutive duplicates exists    
            duplicates++;
            continue;
        }
        pair = new std::pair<double,double>(new_x, new_y);
        (*snapped_curve)[i] = pair;
        i++;
    }

    // padding 
    for(int k=i; k<numPoints; k++) {
        int new_x = PADDING;
        int new_y = PADDING;
        pair = new std::pair<double,double>(new_x, new_y);
        (*snapped_curve)[k] = pair;
    }

    return snapped_curve;
}


std::vector<double> * LSH::getConcatCurve(std::vector<std::pair<double,double>*> *snapped_curve) const {

    std::vector<double> * concat_curve = new std::vector<double>(2 * snapped_curve->size());

    for(int i=0; i< snapped_curve->size(); i++) {
        (*concat_curve)[2*i] = snapped_curve->at(i)->first;
        (*concat_curve)[2*i+1] = snapped_curve->at(i)->second;
        delete snapped_curve->at(i);
    }
    delete snapped_curve;
    return concat_curve;
    //uint64_t key = this->hashFun->amplified_hash(concat_curve, 0);      // get hashed key from amplified function
    //return key;
}

std::vector<double> * LSH::filter(Curve &curve, const double epsilon) const {

    int numPoints = curve.getVector()->size();

    if(numPoints == 1 || numPoints == 2)
        return curve.getVector();

    std::vector<double> *filtered = new std::vector<double>(*(curve.getVector()));
    int i = 0;
    int deleted_points = 0;

    while((i+2) < filtered->size()) {
        double a = filtered->at(i);
        double b = filtered->at(i+1);
        double c = filtered->at(i+2);
        double distance1 = euclideanDistance1D(a,b);
        double distance2 = euclideanDistance1D(b,c);
        // filter out b-point if dist<epsilon
        if(std::max(distance1,distance2) < epsilon) {
            auto it = filtered->begin()+(i+1);
            filtered->erase(it);
            deleted_points++;
        }
        else {
            i++;
        }
    }

    // update the curve with the filtered 
    curve.deleteVector();
    curve.setVector(new std::vector<double>(*filtered));
    
    // padding 
    for(int i=0; i<deleted_points; i++) {
        filtered->push_back(PADDING);
        //(*filtered)[i] = PADDING;
    }

    return filtered;
}

std::vector<double> * LSH::snap_continuous(std::vector<double> *filtered, const double delta, double t) const {

    int numPoints = filtered->size();
    std::pair<double,double>* pair;

    int duplicates = 0, points_inserted = 0;

    for(int i=0; i<numPoints; i++) {
        double y = filtered->at(i);
        double new_y = ((y-t)/delta + 0.5)*delta + t;

        if(i>0 && filtered->at(i-1) == new_y) {   // check if consecutive duplicates exists    
            duplicates++;
            continue;
        }
        else
            points_inserted++;

        (*filtered)[i-duplicates] = new_y;
    }

    // padding 
    for(int i=points_inserted; i<numPoints; i++) {
        (*filtered)[i] = PADDING;
    }

    return filtered;
}

std::vector<double> * LSH::minima_maxima(std::vector<double> *filtered) const {

    int numPoints = filtered->size();

    if(numPoints == 1 || numPoints == 2)
        return filtered;

    int i = 0;
    int deleted_points = 0;

    while((i+2) < filtered->size()) {
        double a = filtered->at(i);
        double b = filtered->at(i+1);
        double c = filtered->at(i+2);
        double lb = std::min(a,c);
        double ub = std::max(a,c);
        // delete b-point if lower_bound < b < upper_bound
        if( (b < ub) && (b > lb) ) {
            auto it = filtered->begin()+(i+1);
            filtered->erase(it);
            deleted_points++;
        }
        else {
            i++;
        }
    }
    
    // padding 
    for(int i=0; i<deleted_points; i++) {
        filtered->push_back(PADDING);
        //(*filtered)[i] = PADDING;
    }

    return filtered;
}

fred::Curve * convertCurve(std::vector<double> *filtered) {
    fred::Curve * curve = new fred::Curve(filtered->size(),1,"");
    for(int i=0; i<filtered->size(); i++) {
        curve->operator[](i).set(0,filtered->at(i));
    }
    return curve;
} 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  GRID  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Grid::Grid(std::vector<Curve *> * CurvesVector, userInput *ui,const double delta) :
    delta(delta),
    t(new std::pair<double,double>(create_rand_double_uniform(delta),create_rand_double_uniform(delta))),
    lsh(new LSH(CurvesVector, ui, delta, t, TYPE_DISCRETE))
{}

Grid::Grid(std::vector<Curve *> * CurvesVector, userInput *ui,const double delta, const double epsilon) :
    delta(delta),
    epsilon(epsilon),
    t(new std::pair<double,double>(create_rand_double_uniform(delta),0)),
    lsh(new LSH(CurvesVector, ui, delta, t->first, epsilon, TYPE_CONTINUOUS))
{}

Grid::~Grid() {
    delete t;
    delete lsh;
}

std::multimap<double, std::pair<Curve *, std::chrono::duration<double>>> & Grid::GridQueryKNN(Curve *CurvePtr, const uint32_t K, // 
      std::multimap<double, std::pair<Curve *, std::chrono::duration<double>> >&map, std::unordered_set<std::string> &explored_set) 
{
    return lsh->queryKNN(CurvePtr, 1, map,explored_set, delta, t, epsilon);
}

bool Grid::rangedSearchCentroidGrid(Curve *curve, std::unordered_set<Curve *> *cluster, std::unordered_set<std::string> &explored_set, const double R, const bool addBarrier, const int typeFrechet) 
{
    return lsh->rangedSearchCentroid(curve, cluster, R, addBarrier, explored_set, delta, *t);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  LSH DISCRETE FRECHET  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


LSH_discrete_frechet::LSH_discrete_frechet(std::vector<Curve *> * CurvesVector, userInput *ui, const double delta) :
    ui(ui)
{
    for(int i=0; i<ui->L; i++) {
        grids.push_back(new Grid(CurvesVector, ui, delta));
    }
}

LSH_discrete_frechet::~LSH_discrete_frechet() {
    for(int i=0; i<grids.size(); i++ )
        delete grids[i];
}


bool LSH_discrete_frechet::queryFile(std::string queryFile) const {

    std::ofstream outputFile(ui->outputFile);
    std::vector<Curve *> *queryCurves;
    queryCurves = parse_file(queryFile);
    if(queryCurves == nullptr)      // if invalid query file is given return false
        return false;
    queryCurves = add_Xaxis(queryCurves);
    std::unordered_set<std::string> explored_set;
    std::multimap<double, std::pair<Point *, std::chrono::duration<double>>> map;
    std::multimap<double, std::pair<Point *, std::chrono::duration<double>>>::iterator it;
    std::multimap<double, std::pair<Point *, std::chrono::duration<double>>>::iterator bruteForceIterator;
    double MAF = 0.0001;
    for(int i=0; i<queryCurves->size(); i++) {
        Point *point = (*queryCurves)[i];
        outputFile << "Query: "+point->getItemId() << '\n'; //<< "Nearest Neighbor-1: "+ std::get<0>(NN) << "\nDistance LSH: "+  std::to_string(std::get<1>(NN))<<"\n\n";
        for(int j=0; j<ui->L; j++)
            map = grids[j]->GridQueryKNN(point, 1, map,explored_set);
        std::multimap<double, std::pair<Point *, std::chrono::duration<double>>> bruteForceMap;
        bruteForceMap = grids[0]->getLSH()->getHTarray()[0]->bruteForce(point, 1, bruteForceMap, TYPE_DISCRETE);
        bruteForceIterator = bruteForceMap.begin();
        int count=0;
        for(it = map.begin(); it != map.end(); ++it) {
            MAF = (it->first/bruteForceIterator->first  > MAF) ? it->first/bruteForceIterator->first : MAF;
            outputFile << "Approximate Nearest neighbor-"<< ++count << ": " << it->second.first->getItemId()<< "\nDistanceLSH: "<<  it->first << "\nDistanceTrue: " << bruteForceIterator->first <<'\n' ;
            outputFile << "tLSH: "<<  it->second.second.count() << "\ntTrue: " << bruteForceIterator->second.second.count() <<'\n' ;
            bruteForceIterator++;
        }
        explored_set.clear();
        map.clear();
        delete point;
    }
    delete queryCurves;
    outputFile << "MAF: "<< MAF << '\n';
    outputFile.close();
    return true;
}

int LSH_discrete_frechet::rangedSearchCentroidFrechet(Curve *curve, std::unordered_set<Curve *> *cluster, const double R, const bool addBarrier) {
    int curves_inserted = 0;
    std::unordered_set<std::string> explored_set;
    for(int i=0; i<grids.size(); i++) {
        curves_inserted += grids[i]->rangedSearchCentroidGrid(curve, cluster, explored_set, R, addBarrier, TYPE_DISCRETE);
    }
    return curves_inserted;                                        
}

void LSH_discrete_frechet::unassignedPoints(std::unordered_map<Point *, std::unordered_set<Point *> *> *clusters, std::unordered_set<Point *> *centroids) {
    for(int i=0; i<grids.size(); i++) 
        grids[i]->getLSH()->unassignedPoints(clusters,centroids);

}





//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  LSH CONTINUOUS FRECHET  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


LSH_continuous_frechet::LSH_continuous_frechet(std::vector<Curve *> * CurvesVector, userInput *ui, const double delta, const double epsilon) :
    ui(ui),
    grid(new Grid(CurvesVector, ui, delta, epsilon))
{}


bool LSH_continuous_frechet::queryFile(std::string queryFile) const {

    std::ofstream outputFile(ui->outputFile);
    std::vector<Curve *> *queryCurves;
    queryCurves = parse_file(queryFile);
    if(queryCurves == nullptr)      // if invalid query file is given return false
        return false;

    std::unordered_set<std::string> explored_set;
    std::multimap<double, std::pair<Point *, std::chrono::duration<double>>> map;
    std::multimap<double, std::pair<Point *, std::chrono::duration<double>>>::iterator it;
    std::multimap<double, std::pair<Point *, std::chrono::duration<double>>>::iterator bruteForceIterator;

    for(int i=0; i<queryCurves->size(); i++) {
        Point *point = (*queryCurves)[i];
        //std::tuple<std::string, double>  NN = queryNN(point);
        outputFile << "Query: "+point->getItemId() << '\n'; //<< "Nearest Neighbor-1: "+ std::get<0>(NN) << "\nDistance LSH: "+  std::to_string(std::get<1>(NN))<<"\n\n";
        
        map = grid->GridQueryKNN(point, ui->N, map,explored_set);
        std::multimap<double, std::pair<Point *, std::chrono::duration<double>>> bruteForceMap;
        bruteForceMap = grid->getLSH()->getHTarray()[0]->bruteForce(point, ui->N, bruteForceMap, TYPE_CONTINUOUS);
        bruteForceIterator = bruteForceMap.begin();
        int count=0;
        for(it = map.begin(); it != map.end(); ++it) {
            outputFile << "Approximate Nearest neighbor-"<< ++count << ": " << it->second.first->getItemId()<< "\nDistanceLSH: "<<  it->first << "\nDistanceTrue: " << bruteForceIterator->first <<'\n' ;
            outputFile << "tLSH: "<<  it->second.second.count() << "\ntTrue: " << bruteForceIterator->second.second.count() <<'\n' ;
            bruteForceIterator++;
        }
        explored_set.clear();
        map.clear();
        delete point;
    }
    delete queryCurves;
    outputFile.close();
    return true;
}

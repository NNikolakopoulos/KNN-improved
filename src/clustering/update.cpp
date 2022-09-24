#include "update.hpp"
#include "../LSH/hash.hpp"
#include "../operations.hpp"
#include "initialization.hpp"


bool update(std::unordered_map<Point*, std::unordered_set<Point*>*>* clusters, std::unordered_set<Point*>* centroids)
{
    double maxDist = DBL_MIN;
    centroids->clear();
    for (auto itr = clusters->begin(); itr != clusters->end(); itr++)
    {
        std::vector<double>* center = mean_vectors(itr->second);
        Point* new_centroid = new Point(std::to_string(0), center);
        centroids->insert(new_centroid);

        double currDist = euclideanDistance( *center , *(itr->first->getVector()));
        if( maxDist <  currDist)
            maxDist = currDist;
    }
    if(maxDist < THRESHOLD )
        return false;

    for (auto itr = clusters->begin(); itr != clusters->end(); itr++) {
        delete itr->second;
        if (itr->first->getItemId() == std::to_string(0))
            delete itr->first;
    }
    clusters->clear();  // initialize clusters for next iteration

    return true;
}


bool update_curves(std::unordered_map<Curve*, std::unordered_set<Curve*>*>* clusters, std::unordered_set<Curve*>* centroids)
{
    double maxDist = DBL_MIN;
    centroids->clear();
    for (auto itr = clusters->begin(); itr != clusters->end(); itr++)
    {
        Curve* c = mean_curve(itr->second);
        Curve* new_centroid = new Point(*c);
        centroids->insert(new_centroid);

        double currDist = discrete_frechet( c , itr->first);
        if( maxDist <  currDist)
            maxDist = currDist;
    }

    if(maxDist < THRESHOLD )
        return false;

    for (auto itr = clusters->begin(); itr != clusters->end(); itr++) {
        delete itr->second;
        if (itr->first->getItemId() == std::to_string(0))
            delete itr->first;
    }
    clusters->clear();  // initialize clusters for next iteration

    return true;
}


std::vector<double>* mean_vectors(std::unordered_set<Point*>* cluster)
{
    auto itr = cluster->begin();
    int size = (*itr)->getVector()->size();
    std::vector<double>* coordinates = new std::vector<double>(size, 0.0);
    while (itr != cluster->end())
    {
        for (int i = 0; i < size; i++)
        {
            (*coordinates)[i] += (*((*itr)->getVector()))[i];
        }
        itr++;
    }
    for (int i = 0; i < size; i++)
    {
        (*coordinates)[i] /= cluster->size();
    }
    return coordinates;
}


std::list<std::tuple<double, double, double, double>> opt_traversal(Curve* c1, Curve* c2)
{
    int length = (c1->getVector()->size())/2;
    int width = (c2->getVector()->size())/2;

    double **C = new double *[length];

    for (int i = 0; i < length; i++)
    {
        C[i] = new double[width];
    }

    double result = 0.0;

    std::vector<double> *v1 = new std::vector<double>(2);
    std::vector<double> *v2 = new std::vector<double>(2);
    

    for (int i = 0; i < length; i++)
    {
        (*v1)[0] = c1->getVector()->at(i);
        (*v1)[1] = c1->getVector()->at(i+length);

        for (int j = 0; j < width; j++)
        {
            (*v2)[0] = c2->getVector()->at(j);
            (*v2)[1] = c2->getVector()->at(j+width);
            double dist = euclideanDistance(*v1, *v2);

            if (i == 0 && j == 0)
            {
                C[i][j] = dist;
            }
            else if (j > 0 && i == 0)
            {
                C[i][j] = std::max(C[i][j - 1], dist);
            }
            else if (i > 0 && j == 0)
            {
                C[i][j] = std::max(C[i - 1][j] , dist);
            }
            else
            {
                double min = std::min(std::min(C[i - 1][j], C[i - 1][j - 1]), C[i][j - 1]);
                C[i][j] = std::max(min , dist);
            }
        }
    }
    delete v1;
    delete v2;
    std::list<std::tuple<double, double, double, double>> traversal;
    int index1 = 2*length - 1;
    int index2 = 2*width - 1;
    int i = length - 1;
    int j = width - 1;
    int min_value;
    traversal.push_front(std::make_tuple(c1->getVector()->at(i), c1->getVector()->at(index1), c2->getVector()->at(j), c2->getVector()->at(index2)));
    while ((index1 > length && index2 > width) && (i>0 && j>0))
    {
        min_value = std::min(std::min(C[i-1][j], C[i][j-1]), C[i-1][j-1]);
        if (min_value == C[i-1][j])
        {
            
            traversal.push_front(std::make_tuple(c1->getVector()->at(--i), c1->getVector()->at(--index1), c2->getVector()->at(j), c2->getVector()->at(index2)));
        }
        else if (min_value == C[i][j-1])
        {
            traversal.push_front(std::make_tuple(c1->getVector()->at(i), c1->getVector()->at(index1), c2->getVector()->at(--j), c2->getVector()->at(--index2)));
        }
        else
        {
            traversal.push_front(std::make_tuple(c1->getVector()->at(--i), c1->getVector()->at(--index1), c2->getVector()->at(--j), c2->getVector()->at(--index2)));
        }
    }
    for (int i = 0; i < length; i++)
    {
        delete[] C[i];
    }

    delete[] C;
    return traversal;
}


Curve* mean_curve(Curve* c1, Curve* c2)
{
    if(c1 == nullptr)
        return c2;
    else if(c2 == nullptr)
        return c1;
    std::list<std::tuple<double, double, double, double>> traversal = opt_traversal(c1, c2);
    int complexity = (c1->getVector()->size())/2;
    int numPointsMeanCurve = traversal.size();
    std::vector<double>* vec = new std::vector<double>(c1->getVector()->size());
    auto itr = traversal.begin();
    double y = (std::get<0>(*itr) + std::get<2>(*itr))/2;
    double x = (std::get<1>(*itr) + std::get<3>(*itr))/2;
    itr++;
    (*vec)[0] = y; 
    (*vec)[complexity] = x; 
    
    for(int i = 1; i<complexity-1; i++) {
        y = (std::get<0>(*itr) + std::get<2>(*itr))/2;
        x = (std::get<1>(*itr) + std::get<3>(*itr))/2;
        itr++;

        if(numPointsMeanCurve > complexity) {
            if(itr == traversal.end()) 
                break;
            double new_y = (std::get<0>(*itr) + std::get<2>(*itr))/2;       //compute next y value of optimal traversal
            double new_x = (std::get<1>(*itr) + std::get<3>(*itr))/2;
            itr++;
            x = (new_x + x)/2;      // get average of two consecutive points in optimal traversal
            y = (new_y + y)/2;
            numPointsMeanCurve--;           // lower the complexity of mean curves, until its equal to the complexity of input curves
        }
        (*vec)[i] = y;   
        (*vec)[i+complexity] = x;  
    }
    Curve* mdfc = new Curve(std::to_string(0), vec);
    return mdfc;
}


Curve *cbt_traversal(CBT *cbt, int index);

Curve* mean_curve(std::unordered_set<Curve*>* cluster)
{
    CBT *cbt = new CBT(*cluster, cluster->size());
    int root = 0;
    return cbt_traversal(cbt, root);
}

Curve *cbt_traversal(CBT *cbt, int index) {
    if(cbt->isLeaf(index))
        return cbt->getCurve(index);
    else {
        Curve *leftCurve = cbt_traversal(cbt, cbt->getLeft(index));
        Curve *RightCurve;
        if(cbt->getCurve(cbt->getRight(index)) != nullptr ) {
            RightCurve =  cbt_traversal(cbt, cbt->getRight(index));
        }
        else
            RightCurve = nullptr;
        
        return mean_curve(leftCurve, RightCurve);
    }
}



/* ~~~~~~~~~~~~~~~       CBT        ~~~~~~~~~~~~~~~~*/

CBT::CBT(std::unordered_set<Curve *> &curves, int numLeaves) {
    int height = std::ceil(std::log2(numLeaves))+1;
    leaf_level_index = pow2(height)/2 - 1;  // index to where the leaf start in array
    size = leaf_level_index + numLeaves;    // size of array

    this->array = new std::vector<Curve *>(size, nullptr);
    auto it = curves.begin();
    for(int i = 0; i< numLeaves; i++) {
        this->array->at(i + leaf_level_index) = *it;
        it++;
        if(it == curves.end())
            break;
    }
}
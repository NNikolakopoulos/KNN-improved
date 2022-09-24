
#include "operations.hpp"


uint32_t pow2(const uint32_t num) {
    uint32_t power = 1;
    for(int i=0; i<num; ++i){
        power *= 2;
    }
    return power;
}

uint64_t mod(const uint64_t a, const uint64_t b) 
{
    uint64_t r = a % b;
    return (r < 0) ? (r + b) : r; 
}

double euclideanDistance(std::vector<double> &a, std::vector<double> &b) 
{
    double dist = 0;
    double temp;
    for(int i=0; i< a.size(); i++) {
        temp = a[i]-b[i];
        dist += temp*temp;
    }
    return sqrt(dist);
}


double euclideanDistance1D(const double a,const double b) 
{
    if(a >= b)
        return a-b;
    else
        return b-a;
}

uint32_t hammingDistance(const uint64_t n1, const uint64_t n2)
{
    int x = n1 ^ n2;
    uint32_t nbits = 0;
 
    while (x > 0) {
        nbits = nbits + (x & 1) ;
        x = x >> 1;
    }
 
    return nbits;
}

std::queue<uint32_t> & hamDist1(std::queue<uint32_t> & queue, const uint32_t key, std::unordered_set<uint32_t> & expl_set_hamDist, const uint32_t length) {

    std::string binary_str = std::bitset<32>(key).to_string(), tempStr;
    uint32_t newKey;

    binary_str = binary_str.substr(32-length,length);
    for(int i = 0; i < length; ++i) {

        tempStr = binary_str;

        if(binary_str.at(i) == '1') 
            tempStr.replace(i,1,"0");
        else if(binary_str.at(i) == '0')
            tempStr.replace(i,1,"1");

        newKey = stol(tempStr,nullptr,2);                                       // convert binary string to int

        if(expl_set_hamDist.find(key) != expl_set_hamDist.end())          // if this numbered is already explored then dont insert into the queue
            continue;
        else {
            queue.push(newKey);
            expl_set_hamDist.insert(newKey);
        }
    }
    return queue;
}

//returns the index of the first number in the vector that is bigger than "num"
int search_num_range(std::vector<double> &vec, double &num)
{
    if (vec.size() == 1)
    {
        return vec[0];
    }
    int low = 0;
    int high = vec.size()-1;
    int temp = (high+low)/2;
    while (true)
    {
        if (low == high - 1)
        {
            return high;
        }
        if (num > vec[temp])
        {
            low = temp;
            temp = (low+high)/2;
        }
        else
        {
            high = temp;
            temp = (low+high)/2;
        }
    }
}

double compute_delta() {
    //double smallNum = 10e-3;
    //return 4 * curve->getVector()->size() * smallNum;
    return 0.5;
}
double discrete_frechet(Curve *c1, Curve *c2) {

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

    result = C[length - 1][width - 1];

    for (int i = 0; i < length; i++)
    {
        delete[] C[i];
    }

    delete[] C;
    delete v1;
    delete v2;

    return result;
} 

//the below code  is not used
double dynamic_tw(Curve *c1, Curve *c2) {

    int length = (c1->getVector()->size())/2;
    int width = (c2->getVector()->size())/2;

    double **C = new double *[length];

    for (int i = 0; i < length; i++)
    {
        C[i] = new double[width];
    }

    double minimum, result = 0.0;

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

            if (i == 0 && j == 0)
            {
                C[i][j] = euclideanDistance(*v1, *v2);
            }
            else if (j > 0 && i == 0)
            {
                C[i][j] = C[i][j - 1] + euclideanDistance(*v1, *v2);
            }
            else if (i > 0 && j == 0)
            {
                C[i][j] = C[i - 1][j] + euclideanDistance(*v1, *v2);
            }
            else
            {
                minimum = std::min(std::min(C[i - 1][j], C[i - 1][j - 1]), C[i][j - 1]);
                C[i][j] = minimum + euclideanDistance(*v1, *v2);
                minimum = 0.0;
            }
        }
    }

    result = C[length - 1][width - 1];

    for (int i = 0; i < length; i++)
    {
        delete[] C[i];
    }

    delete[] C;
    delete v1;
    delete v2;

    return result;
}


double traverse(double **C, Curve*c1, Curve *c2,int i, int j,std::vector<double> *v1,std::vector<double> *v2);

double optimal_traversal(Curve *c1, Curve *c2) {

    int length = (c1->getVector()->size())/2;
    int width = (c2->getVector()->size())/2;
    double minimum, result = 0.0;


    //std::vector<double> *traversal = new std::vector<double>(std::max(length,width));

    double **C = new double *[length];

    for (int i = 0; i < length; i++)
    {
        C[i] = new double[width];
    }
    //std::memset(C, -1.0, (length*width)*sizeof(C[0][0]));
    for(int i =0; i<length; i++)
        for(int j=0; j<width; j++) 
            C[i][j] = -1;

    std::vector<double> *v1 = new std::vector<double>(2);
    std::vector<double> *v2 = new std::vector<double>(2);

    int i = length-1;
    int j = width-1;

    
    int result1 = traverse(C, c1,c2, i,  j,v1,v2);
    result = C[length - 1][width - 1];


    for (int i = 0; i < length; i++)
    {
        delete[] C[i];
    }

    delete[] C;
    delete v1;
    delete v2;

    return result;
}

double traverse(double **C, Curve *c1,Curve *c2,int i, int j,std::vector<double> *v1,std::vector<double> *v2) {


    // if its already computed , then return this value
    if(C[i][j] != -1)
        return C[i][j];

    //else we must compute 
    
    (*v1)[0] = c1->getVector()->at(i);
    (*v1)[1] = c1->getVector()->at(i+((c1->getVector()->size())/2));
    (*v2)[0] = c2->getVector()->at(j);
    (*v2)[1] = c2->getVector()->at(j+((c2->getVector()->size())/2));
    if(i == 0 && j == 0) {
        C[i][j] = euclideanDistance(*v1, *v2);
    }
    else if (j > 0 && i == 0)
    {
        C[i][j] = traverse(C,c1,c2,i,j-1,v1,v2) + euclideanDistance(*v1, *v2);
    }
    else if (i > 0 && j == 0)
    {
        C[i][j] = traverse(C,c1,c2,i-1,j,v1,v2) + euclideanDistance(*v1, *v2);
    }
    else
    {
        int minimum = std::min(std::min(traverse(C,c1,c2,i-1,j,v1,v2),traverse(C,c1,c2,i-1,j-1,v1,v2)), traverse(C,c1,c2,i,j-1,v1,v2));
        //int minimum = std::min(std::min(C[i - 1][j], C[i - 1][j - 1]), C[i][j - 1]);
        C[i][j] = minimum + euclideanDistance(*v1, *v2);
        minimum = 0.0;
    }
    return C[i][j];
}
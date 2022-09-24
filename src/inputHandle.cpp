#include "inputHandle.hpp"
int get_user_input_commant_promt(userInput *& usrIn, const int argc, char **argv) {

    argv[argc-1][strlen(argv[argc-1])] = '\0';             // last char is '\n'

    std::string tempStr;
    usrIn = new userInput;
    int filesCounter = 0;

    usrIn->k = 4;
    usrIn->L = 5;
    usrIn->N = 1;
    usrIn->R = 10000;
    usrIn->M = 10;
    usrIn->probes = 2;

    for(int i=1; i<argc; i++) {

        if( argv[i][0] == '-') {

            // if there is nothing after the flag, return error
            if( i+1 >= argc)
                return -1;
            // if after the flag is another flag, return error
            else if( argv[i+1][0] == '-') 
                return -1;

            else if(strcmp(argv[i],"-i") == 0) {
                usrIn->inputFile += argv[i+1];
                filesCounter++;
            }
            else if(strcmp(argv[i],"-q") == 0) {
                usrIn->queryFile += argv[i+1];
                filesCounter++;
            }
            else if(strcmp(argv[i],"-o") == 0) {
                usrIn->outputFile += argv[i+1];
                filesCounter++;
            }
            else if(strcmp(argv[i], "-k") == 0)
                usrIn->k = atoi(argv[i+1]);
            else if(strcmp(argv[i], "-L") == 0) 
                usrIn->L = atoi(argv[i+1]);
            else if(strcmp(argv[i],"-N") == 0) 
                usrIn->N = atoi(argv[i+1]);
            else if(strcmp(argv[i], "-M") == 0)
                usrIn->M = atoi(argv[i+1]);
            else if(strcmp(argv[i],"-probes") == 0)
                usrIn->probes = atoi(argv[i+1]);
            else if(strcmp(argv[i],"-R") == 0)
                usrIn->R = atoi(argv[i+1]);
            else if(strcmp(argv[i],"-algorithm") == 0)
                usrIn->algorithm = argv[i+1];
            else if(strcmp(argv[i],"-metric") == 0)
                usrIn->metric = argv[i+1];
            else if(strcmp(argv[i],"-delta") == 0)
                usrIn->delta = atof(argv[i+1]);
        }
    }
    if(filesCounter < 3) {
        std::cout << "Please give all necessary file names\n";
        return -1;
    }           
    return 1; // all good
} 


int user_input_cluster(clusterInput *& usrIn, const int argc, char **argv) {

    argv[argc-1][strlen(argv[argc-1])] = '\0';             // last char is '\n'

    std::string tempStr;
    usrIn = new clusterInput;
    int filesCounter = 0;
    usrIn->silhouette = false;

    for(int i=1; i<argc; i++) {

        if( argv[i][0] == '-') {

            if(strcmp(argv[i],"-i") == 0) {
                usrIn->inputFile += argv[i+1];
                filesCounter++;
            }
            else if(strcmp(argv[i],"-c") == 0) {
                usrIn->configFile += argv[i+1];
                filesCounter++;
            }
            else if(strcmp(argv[i],"-o") == 0) {
                usrIn->outputFile += argv[i+1];
                filesCounter++;
            }
            else if(strcmp(argv[i], "-complete") == 0) {
                usrIn->complete = true;
            }
            else if(strcmp(argv[i], "-update") == 0) {
                usrIn->update = argv[i+1];
            }
            else if(strcmp(argv[i], "-assignment") == 0) {
                usrIn->assignment = argv[i+1];
            }
            else if(strcmp(argv[i], "-silhouette") == 0) {
                usrIn->silhouette = true;
            }
        }
    }
    if(filesCounter < 3) {
        std::cout << "Please give all necessary file names\n";
        return -1;
    }           
    return 1; // all good
} 

std::vector<Point *> * parse_file(std::string fileName) {

    Point *point;
    std::vector<double> *pointVec;
    std::vector<Point *> *pointsVec = new std::vector<Point *>;
    std::ifstream inFileStream(fileName);
    std::string line, token, itemId;
    
    while(getline(inFileStream, line, '\r')) {
        if(line == "\n")                            // thats the last empty line
            break;

        std::string token;
        std::stringstream s(line);
        getline(s, itemId, '\t');                     // get the item-id of the vector 
        if(itemId[0]=='\n') itemId.erase(0,1);
        pointVec = new std::vector<double>;
        while(getline(s, token, '\t')) {              // now get the coordinates
            pointVec->push_back(std::stod(token));
        }
        point = new Point(itemId, pointVec);
        pointsVec->push_back(point);               // save this point to insert later into Hash Table
        //pointVec.clear();
    }

    uint32_t tableSize = pointsVec->size();
    if(tableSize == 0) {
        std::cout << "Please give a valid file!\n";
        delete pointsVec;
        return NULL;
    }
    inFileStream.close();
    return pointsVec;
}
   
//IN CASE OF BAD/INVALID CONFIG FILE, THIS FUNCTION RETURNS nullptr
std::vector<uint32_t> * getConfigFile(std::string fileName) {

    std::vector<uint32_t> *vec = new std::vector<uint32_t>(6);
    (*vec)[0] = 0; (*vec)[1] = 3; (*vec)[2] = 4; (*vec)[3] = 10; (*vec)[4] = 3; (*vec)[5] = 2; 


    std::string line;
    std::ifstream inFileStream(fileName);
    while(getline(inFileStream, line, '\n')) {
        if(line == "\n")                            // thats the last empty line
            break;
        std::string token, previousToken;
        std::stringstream s(line);
        if(getline(s, token, ' ')) {
            previousToken = token;
            if(getline(s, token, '\n')) {
                if(previousToken == "number_of_clusters:") 
                    (*vec)[0] = std::stoi(token);
                else if(previousToken == "number_of_vector_hash_tables:")
                    (*vec)[1] = std::stoi(token);
                else if(previousToken == "number_of_vector_hash_functions:")
                    (*vec)[2] = std::stoi(token);
                else if(previousToken == "max_number_M_hypercube:")
                    (*vec)[3] = std::stoi(token);
                else if(previousToken == "number_of_hypercube_dimensions:")
                    (*vec)[4] = std::stoi(token);
                else if(previousToken == "number_of_probes:")
                    (*vec)[5] = std::stoi(token);

            }
            else if(previousToken == "number_of_clusters:")
                return nullptr;    
        }
        else
            return nullptr;
    }

    inFileStream.close();
    return vec;
}
void output(clusterInput* usrIn, std::unordered_map<Point*, std::unordered_set<Point*>*>* clusters, std::chrono::duration<double> duration, std::vector<double> s)
{
    std::ofstream outputFile(usrIn->outputFile);
    outputFile << "Algorithm: " << usrIn->assignment << " with " << usrIn->update << std::endl;

    auto itr = clusters->begin();
    for (int i = 0; i < clusters->size(); i++)
    {
        outputFile << "CLUSTER-" << i+1 << " {size: " << itr->second->size() << ", centroid: ";
        std::vector<double>* v = itr->first->getVector();
        for (double n : *v)
        {
            outputFile << n << ' ';
        }
        outputFile << '}' << std::endl;
        itr++;
    }
    outputFile << "clustering_time: " << duration.count() << " seconds" << std::endl;
    if (usrIn->silhouette == true)
    {
        outputFile << "Silhouette: [";
        for (double n : s)
        {
            outputFile << n << ' ';
        }
        outputFile << ']' << std::endl;
    }
    itr = clusters->begin();
    if (usrIn->complete == true)
    {
        outputFile << '\n' << std::endl;
        for (int i = 0; i < clusters->size(); i++)
        {
            outputFile << "CLUSTER-" << i+1 << " {Points: ";
            for (auto itr2 = itr->second->begin(); itr2 != itr->second->end(); itr2++)
            {
                outputFile << (*itr2)->getItemId() << ' ';
            }

            outputFile << '}' << std::endl;
            itr++;
        }
    }
    
    outputFile.close();
}

std::vector<Curve *> * add_Xaxis(std::vector<Curve *> *CurvesVector) {
    int numPoints = CurvesVector->at(0)->getVector()->size();           // get the number of points in each vector
    for(int j=0; j<CurvesVector->size(); j++)
        for(int i=0; i<numPoints; i++)
            (*CurvesVector)[j]->add(i);
    return CurvesVector;    
}
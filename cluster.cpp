#include "headers.hpp"
#include "src/inputHandle.hpp"
#include "src/LSH/hash.hpp"
#include "src/clustering/initialization.hpp"
#include "src/clustering/assignment.hpp"
#include "src/clustering/reverse.hpp"
#include "src/clustering/update.hpp"


int main(int argc, char **argv)
{
    clusterInput *usrIn;
    if(user_input_cluster(usrIn,argc,argv) == -1 ) {
        std::cout << "Please give valid arguments!\n";
    }
    std::vector<Point *> * PointsVector = parse_file(usrIn->inputFile);
    if(PointsVector == NULL) return -1;
    std::vector<uint32_t> * parameters = getConfigFile(usrIn->configFile);
    std::unordered_set<Point*>* centroids;
    std::unordered_map<Point*, std::unordered_set<Point*>*>* clusters = new std::unordered_map<Point*, std::unordered_set<Point*>*>;

    userInput *ui = new userInput;
    int flagFrechet = 0;
    ui->k = (*parameters)[2];
    ui->inputFile = usrIn->inputFile;
    ui->outputFile = usrIn->outputFile;
    ui->L = (*parameters)[1];
    ui->M = (*parameters)[3];
    ui->probes = (*parameters)[5];

    auto start = std::chrono::high_resolution_clock::now();
    reverseAssignment* revAs = nullptr;
    if (usrIn->update == "Mean_Vector" || usrIn->update == "mean_vector")
    {
        centroids = initialization(PointsVector, (*parameters)[0]);
        if (usrIn->assignment == "LSH" || usrIn->assignment == "lsh")
        {
            revAs = new reverseAssignment(PointsVector, ui, centroids, 1);
        }
        else if (usrIn->assignment == "Hypercube" || usrIn->assignment == "hypercube")
        {
            revAs= new reverseAssignment(PointsVector, ui, centroids, 0);
        }
        do
        {
            assignment(PointsVector, centroids, clusters, usrIn->assignment, ui, revAs);
        }
        while(update(clusters, centroids));
    }
    else if (usrIn->update == "Mean_Frechet" || usrIn->update == "mean_frechet")
    {
        flagFrechet = TYPE_DISCRETE;
        PointsVector = add_Xaxis(PointsVector);
        centroids = initialization_curves(PointsVector, (*parameters)[0]);
        if (usrIn->assignment == "LSH_Frechet" || usrIn->assignment == "lsh_frechet")
        {
            revAs = new reverseAssignment(PointsVector, ui, centroids, 1);
        }
        do
        {
            assignment_curves(PointsVector, centroids, clusters, usrIn->assignment, ui, revAs);
        }
        while(update_curves(clusters, centroids));
    }
    else {
        return -1;
    }
    
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = stop - start;
    std::vector<double> s;
    flagFrechet = 0;
    if (usrIn->silhouette == true)
    {
        s = silhouette(clusters, flagFrechet);
    }
    output(usrIn, clusters, duration, s);

    deleteClusters(clusters);
    deletePointsVector(PointsVector, usrIn->assignment);
    delete revAs;
    delete centroids;
    delete usrIn;
    delete parameters;
    return 0;
}
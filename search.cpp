#include "headers.hpp"
#include "./src/inputHandle.hpp"
#include "./src/LSH/hash.hpp"
#include "./src/hypercube/hyperc.hpp"

int main(int argc, char **argv) {

    // get user input
    userInput *usrIn;
    if(get_user_input_commant_promt(usrIn,argc,argv) == -1 ) {
        std::cout << "Please give valid arguments!\n";
        return -1;
    }

    std::vector<Curve *> * CurvesVector = parse_file(usrIn->inputFile);
    if(CurvesVector == NULL) return -1;
    std::string input;
    double epsilon = 0.01;
    if (usrIn->algorithm == "LSH" || usrIn->algorithm == "lsh")
    {
        LSH l(CurvesVector, usrIn, TYPE_VECTOR);
        while( true ) {
            if( !(l.queryFile(usrIn->queryFile)) ) {
                std::cin >> usrIn->queryFile;
                continue;
            }

            std::cout << "\nDo you want to repeat the process for another query file? ( y / n )\n";
            std::cin >> input;

            //stop the loop
            if(input == "n") 
                break;

            std::cout << "Please insert name of new file:\n";
            std::cin >> input;
            usrIn->queryFile = input;
        }
    }
    else if (usrIn->algorithm == "Frechet" || usrIn->algorithm == "frechet")
    {
        CurvesVector = add_Xaxis(CurvesVector);
        if (usrIn->metric == "discrete")
        {
            LSH_discrete_frechet l(CurvesVector, usrIn, usrIn->delta);
            while( true ) {
                if( !(l.queryFile(usrIn->queryFile)) ) {
                    std::cin >> usrIn->queryFile;
                    continue;
                }

                std::cout << "\nDo you want to repeat the process for another query file? ( y / n )\n";
                std::cin >> input;

                //stop the loop
                if(input == "n") 
                    break;

                std::cout << "Please insert name of new file:\n";
                std::cin >> input;
                usrIn->queryFile = input;
            }
        }
        else if (usrIn->metric == "continuous")
        {
            LSH_continuous_frechet l(CurvesVector, usrIn, usrIn->delta, epsilon);
            while( true ) {
                if( !(l.queryFile(usrIn->queryFile)) ) {
                    std::cin >> usrIn->queryFile;
                    continue;
                }

                std::cout << "\nDo you want to repeat the process for another query file? ( y / n )\n";
                std::cin >> input;

                //stop the loop
                if(input == "n") 
                    break;

                std::cout << "Please insert name of new file:\n";
                std::cin >> input;
                usrIn->queryFile = input;
            }
        }
        else
        {
            std::cout << "Wrong metric given!" << std::endl;
            return -1;
        }
    }
    else if (usrIn->algorithm == "Hypercube" || usrIn->algorithm == "hypercube")
    {
        RandomizedProjection rp(CurvesVector, usrIn);
        while( true ) {
            
            if( !(rp.queryFile(usrIn->queryFile)) ) {
                std::cin >> usrIn->queryFile;
                continue;
            }
            std::cout << "\nDo you want to repeat the process for another query file? ( y / n )\n";
            std::cin >> input;

            //stop the loop
            if(input == "n") 
                break;

            std::cout << "Please insert name of new file:\n";
            std::cin >> input;
            usrIn->queryFile = input;
        }
    }
    std::cout << "Exiting...\n";
    delete CurvesVector;
}
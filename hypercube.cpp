#include "headers.hpp"
#include "./src/inputHandle.hpp"
#include "./src/hypercube/hyperc.hpp"

int main(int argc, char **argv) {

    // get user input
    userInput *usrIn;
    if(get_user_input_commant_promt(usrIn,argc,argv) == -1 ) {
        std::cout << "Please give valid arguments!\n";
    }
    std::vector<Point *> * PointsVector = parse_file(usrIn->inputFile);
    if(PointsVector == NULL) return -1;
    
    std::string input;
    RandomizedProjection rp(PointsVector, usrIn);
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

    std::cout << "Exiting...\n";
    delete PointsVector;
}


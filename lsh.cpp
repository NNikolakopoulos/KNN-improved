#include "headers.hpp"
#include "./src/inputHandle.hpp"
#include "./src/LSH/hash.hpp"

int main(int argc, char **argv) {

    // get user input
    userInput *usrIn;
    if(get_user_input_commant_promt(usrIn,argc,argv) == -1 ) {
        std::cout << "Please give valid arguments!\n";
        return -1;
    }

    std::vector<Curve *> * CurvesVector = parse_file(usrIn->inputFile);
    CurvesVector = add_Xaxis(CurvesVector);
    if(CurvesVector == NULL) return -1;
    std::string input;
    double delta = compute_delta();
    double epsilon = 0.1;
    //LSH l(CurvesVector, usrIn, TYPE_VECTOR);
    //LSH_discrete_frechet l(CurvesVector, usrIn, delta);
    LSH_continuous_frechet l(CurvesVector, usrIn, delta, epsilon);
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

    std::cout << "Exiting...\n";
    delete usrIn;
    for(int i=0; i<CurvesVector->size();i++) {
        delete (*CurvesVector)[i];
    }
    delete CurvesVector;

}


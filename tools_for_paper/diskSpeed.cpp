#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "tclap/CmdLine.h" // command line arguments

using namespace std;
using namespace TCLAP;     // command line arguments

int main(int argc, char * argv[]) {
    
    //processing command line options
    int numberOfSequences, lengthOfSequence;
    string inputFile;
    try{
        CmdLine cmd("trim with z allowed low quality nucleotides per row", ' ', "1.0", true);
        ValueArg<int>    rowsArg(   "r", "reads",  "number of reads",     true,  0,  "integer", cmd);
        ValueArg<int>    lengthArg( "l", "length", "length of each read", true,  0,  "integer", cmd);
        ValueArg<string> infileArg ("i", "infile", "input file name",     true,  "", "string",  cmd);
        cmd.parse( argc, argv );
        numberOfSequences = rowsArg.getValue();
        lengthOfSequence  = lengthArg.getValue();
        inputFile         = infileArg.getValue();
        
    } catch (ArgException &e) {
        cerr << "ARGUMENT ERROR: " << e.error() << " for arg " << e.argId() << endl;
        return EXIT_FAILURE;
    }
    
    // open file
    ifstream in(inputFile, ios::in);
    string zeile;
    
    // do something very easy with the file
    unsigned long long int sum = 0;
    int row_sum;
    for (int z = 0; z < numberOfSequences; z++) {
        in.ignore(numeric_limits<streamsize>::max(), '\n');
        in.ignore(numeric_limits<streamsize>::max(), '\n');
        in.ignore(numeric_limits<streamsize>::max(), '\n');
        getline(in,zeile);
        // do something very easy with the row
        row_sum = 0;
        for (int i=0; i<lengthOfSequence; i++) {
            row_sum += zeile[i];
        }
        sum += (row_sum % 10);
    }
    
    cout << sum << endl;
    
    return EXIT_SUCCESS;
    
}

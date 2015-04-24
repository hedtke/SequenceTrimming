#include <iostream>
#include <fstream>
#include <cstdlib>

// Ivo Hedtke
// 24. April 2015

// arg1: number or reads (lines)
// arg2: length of reads (#cols)

using namespace std;

int main(int argc, const char * argv[]) {
    int reads = atoi(argv[1]);
    int length = atoi(argv[2]);
    
    ofstream myfile;
    myfile.open ("random.fastq", ios::out);
    
    for (int row=0; row<reads; row++) {
        myfile << "idline\n";
        myfile << "CGAT\n";
        myfile << "idline\n";
        for (int col=0; col<length; col++) {
            int ascii_index = 44 + (rand() % 30);
            myfile << (char) ascii_index;
        }
        myfile << "\n";
    }
}

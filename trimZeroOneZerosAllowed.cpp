/*******************************************************************************
 *
 * trimZeroOneZerosAllowed.cpp
 *
 * DESCRIPTION: Given an ASCII grid of zeros and ones or a FASTQ file and a threshold.
 *              n lines. l numbers per line. Parameter k. Find two borders
 *              0 <= l1 <= l2 <= l-1, such that the window on the n x l grid starting
 *              at column l1 and ending at column l2 consists only of ones (quality >=
 *              threshold) and at most k zeros per line. Here it is allowed to ignore
 *              some lines. Determine l1, l2 and the selected (not ignored) lines in
 *              such a way, that the area
 *                 (width*height = (l2 - l1 - 1) * (number of selected lines)
 *              is maximized. If the option for the output file is given, a CSV is
 *              written with the number of selected lines for all given pairs of
 *              left and right borders.
 *
 * RUNTIME: O( n*l + l^2 )
 *
 * AUTHORS: Ivo Hedtke (hedtke@informatik.uni-halle.de)
 *          Matthias Mueller-Hannemann (muellerh@informatik.uni-halle.de)
 *
 * CREATED: 21 Feb 2013
 *
 * LAST CHANGE: 16 Dec 2013
 *
 */

#include "tclap/CmdLine.h" // command line arguments
#include "ComputeMatrices.h" // new algorithms
#include "Results.h" // output on screen or in CSV

using namespace std;
using namespace TCLAP; // command line arguments
using namespace ComputeMatrices; // new algorithms
using namespace Results; // output on screen or in CSV

int main(int argc, char * argv[]) {

    try{

        // read command line parameters
        CmdLine cmd("trim with z allowed zeros per row", ' ', "1.0", true);
        ValueArg<int> rowsArg("r", "rows", "number of sequences/rows", true, 0, "integer", cmd);
        ValueArg<int> lengthArg("l", "length", "length of each sequence/row", true, 0, "integer", cmd);
        ValueArg<int> zerosArg("z", "zeros", "number of allowed zeros per sequence/row", true, 0, "integer", cmd);
        ValueArg<string> infileArg ("i", "infile",  "input file name", true,  "", "string", cmd);
        ValueArg<string> outfileArg("o", "outfile", "output file name (CSV format)", false, "", "string", cmd);
        SwitchArg fastqSwitch("f", "fastq", "input file is in FASTQ format", cmd, false);
        ValueArg<int> thresholdArg("t", "threshold", "quality is ok if quality score >= threshold", false, -1, "integer", cmd);
        ValueArg<int> shiftArg("s", "shift", "shift for char -> quality conversion", false, -1, "integer", cmd);
        cmd.parse( argc, argv );
        int numberOfSequences = rowsArg.getValue();
        int lengthOfSequence = lengthArg.getValue();
        int numberOfAllowedZerosPerSequence = zerosArg.getValue();
        string inputFile = infileArg.getValue();
        string outputFile = outfileArg.getValue();
        bool fastq = fastqSwitch.getValue();
        int threshold = thresholdArg.getValue();
        int shift = shiftArg.getValue();

        // FASTQ files need a threshold
        if ( fastq && (!thresholdArg.isSet() || !shiftArg.isSet()) ){
            cerr << "FASTQ files need a threshold and a shift" << endl;
            return EXIT_FAILURE;
        }
        if ( thresholdArg.isSet() != shiftArg.isSet() ) {
            cerr << "threshold and shift can only be used together, one is missing" << endl;
        }

        // compute matrix c_z for z-zeros
        int** c = trimZeroOneZerosAllowed(inputFile,
                                          numberOfSequences,
                                          lengthOfSequence,
                                          numberOfAllowedZerosPerSequence,
                                          fastq,
                                          threshold,
                                          shift);

        // output in CSV or on terminal
        if (outfileArg.isSet()) {
            exportResultMatrixAsCSV(c, lengthOfSequence, outputFile);
        } else {
            printMaxArea(c, lengthOfSequence, numberOfSequences);
        }
    
        free_matrix(c, lengthOfSequence);

    } catch (ArgException &e) {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
    }

    return EXIT_SUCCESS;

}

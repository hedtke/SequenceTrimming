/*******************************************************************************
 *
 * trimZeroOnePercentZerosAllowed.cpp
 *
 * DESCRIPTION: Given a FASTQ file and a threshold. n lines. l numbers per line.
 *              Parameter p. Find two borders 0 <= l1 <= l2 <= l-1, such that
 *              the window on the (n x l) grid starting at column l1 and ending
 *              at column l2 consists only of high quality nucleotides (quality
 *              score >= threshold) and at most p percent low quality
 *              nucleotides per line. Here it is allowed to ignore some lines.
 *              Determine l1, l2 and the selected (not ignored) lines in such a
 *              way, that the area
 *                 (width*height = (l2 - l1 - 1) * (number of selected lines)
 *              is maximized. If the option for the output file is given, a CSV
 *              is written with the number of selected lines for all given pairs
 *              of left and right borders.
 *
 * RUNTIME: O( l^2 * n )
 *
 * AUTHORS: Ivo Hedtke (ivo.hedtke@uni-osnabrueck.de)
 *          Matthias Mueller-Hannemann (muellerh@informatik.uni-halle.de)
 *
 * CREATED: 19 Dec 2013
 *
 * LAST CHANGE: 29 Jul 2014
 *
 */

#include "tclap/CmdLine.h"           // command line arguments
#include "ComputeMatrices.h"         // trimming algorithms
#include "ComputeMatricesParallel.h" // parallel trimming algorithms
#include "Results.h"                 // output on screen or in CSV

using namespace std;
using namespace TCLAP;           // command line arguments
using namespace ComputeMatrices; // trimming algorithms
using namespace Results;         // output on screen or in CSV

int main(int argc, char * argv[]) {
    
    //START: processing command line options
    int numberOfSequences, lengthOfSequence, threshold, shift, numThreads;
    string inputFile, outputFile;
    double percentOfAllowedZerosPerSequence;
    
    try{
        
        // read command line parameters
        CmdLine cmd("trim with p percent allowed low quality nucleotides per row", ' ', "1.0", true);
        ValueArg<int>    rowsArg(      "r", "reads",       "number of reads",                                          true,  0,   "integer", cmd);
        ValueArg<int>    lengthArg(    "l", "length",      "length of each read",                                      true,  0,   "integer", cmd);
        ValueArg<double> percentArg(   "p", "percent",     "percent of allowed zeros per read: value between 0 and 1", true,  0.0, "double",  cmd);
        ValueArg<string> infileArg(    "i", "infile",      "input file name",                                          true,  "",  "string",  cmd);
        ValueArg<string> outfileArg(   "o", "outfile",     "output file name (CSV format)",                            false, "",  "string",  cmd);
        ValueArg<int>    thresholdArg( "t", "threshold",   "quality is ok if quality score >= threshold",              true,  -1,  "integer", cmd);
        ValueArg<int>    shiftArg(     "s", "shift",       "shift for char -> quality conversion",                     true,  -1,  "integer", cmd);
        ValueArg<int>    numThreadsArg("w", "workthreads", "number of parallel worker threads",                        false,  0,  "integer", cmd);
        
        cmd.parse( argc, argv );
        numberOfSequences                = rowsArg.getValue();
        lengthOfSequence                 = lengthArg.getValue();
        percentOfAllowedZerosPerSequence = percentArg.getValue();
        inputFile                        = infileArg.getValue();
        outputFile                       = outfileArg.getValue();
        threshold                        = thresholdArg.getValue();
        shift                            = shiftArg.getValue();
        numThreads                       = numThreadsArg.getValue();
        
    } catch (ArgException &e) {
        cerr << "ERROR: " << e.error() << " for arg " << e.argId() << endl;
        return EXIT_FAILURE;
    }
    //END: processing command line options
    
    //START: now compute optimal trimming parameters
    vector<vector<int>> c; // compute matrix c_p for p-percent
    
    if (numThreads == 0){// sequential mode
        c = trimZeroOnePercentZerosAllowed(inputFile,numberOfSequences,
                                           lengthOfSequence,
                                           percentOfAllowedZerosPerSequence,
                                           threshold,shift);
    } else {// parallel mode
        c = trimZeroOnePercentZerosAllowedPar(inputFile,numberOfSequences,
                                              lengthOfSequence,
                                              percentOfAllowedZerosPerSequence,
                                              threshold,shift,numThreads );
    }
    //END: now compute optimal trimming parameters
    
    //START: output in CSV or on terminal
    if (outputFile != "") {
        exportMatrix(c,outputFile);
    } else {
        printMaxArea(c, numberOfSequences);
    }
    //END: output in CSV or on terminal
    
    return EXIT_SUCCESS;
    
}

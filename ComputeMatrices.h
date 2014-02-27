/*******************************************************************************
 *
 * ComputeMatrices.h
 *
 * DESCRIPTION: Implementation of the algorithms for the problems:
 *              0-zeros:   trimZeroOne                    (line 30)
 *              z-zeros:   trimZeroOneZerosAllowed        (line 103)
 *              p-percent: trimZeroOnePercentZerosAllowed (line 254)
 *              m-mean:    trimIntegerMean                (line 308)
 *
 * RUNTIMES: If the input has r reads of length l:
 *           0-zeros:   O( r * l )
 *           z-zeros:   O( r * l )
 *           p-percent: O( r * l^2 )
 *           m-mean:    O( r * l^2 )
 *
 * AUTHORS: Ivo Hedtke (ivo.hedtke@uni-osnabrueck.de)
 *          Matthias Mueller-Hannemann (muellerh@informatik.uni-halle.de)
 *
 * CREATED: 13 Dec 2013
 *
 * LAST CHANGE: 27 Feb 2014
 *
 */

#include <fstream>
#include <string>
#include <vector>

using namespace std;

namespace ComputeMatrices {

    // 0-zeros
    // =======
    // readFASTQ => each time skip 3 lines
    // (thresholdGoodValues >= 0) => lines of the grid are not '0' and '1' so use a threshold:
    //                               (quality < threshold)? '0' : '1'
    vector<vector<int> > trimZeroOne(
        const string& inputfile,
        const int& numberOfSequences,
        const int& lengthOfSequence,
        const int& thresholdGoodValues,
        const int& shiftToConvertChars)
    {

        int thresholdPlusShift = thresholdGoodValues + shiftToConvertChars;

        // alloc and init c and cT
        // cT = counter of triangles
        // c(i,j)=m means, there are m lines where a 1-block starts at i & ends at j
        vector<vector<int> > c (lengthOfSequence, vector<int>(lengthOfSequence,0));
        vector<vector<int> > cT (lengthOfSequence, vector<int>(lengthOfSequence,0));

        // read the file row by row
        string zeile;

        // open file
        ifstream in(inputfile, ios::in);

        bool stillInOneBlock = false;
        int startOfOneBlock;
        for (int z = 0; z < numberOfSequences; z++) {
            getline(in,zeile); // skip 3 lines
            getline(in,zeile); // skip 3 lines
            getline(in,zeile); // skip 3 lines
            getline(in,zeile);
            startOfOneBlock = 0;
            zeile.append("!"); // dummy 0 at the end
            // ASSERT: ASCII("!")=33 and this is the smallest possible char in a quality score string
            // we need a dummy "bad" quality score at the end for our algorithm
            for (int i = 0; i <= lengthOfSequence; i++) {
                if (zeile.at(i) < thresholdPlusShift) {
                    if (stillInOneBlock) {
                        stillInOneBlock = false;
                        cT.at(startOfOneBlock).at(i-1)++;
                    }
                } else {
                    if (!stillInOneBlock) {
                        stillInOneBlock = true;
                        startOfOneBlock = i;
                    }
                }
            }
        }

        // compute c from cT
        vector<int> columnSumAbove (lengthOfSequence,0);
        // first fill the last column of c
        c.at(0).at(lengthOfSequence-1) = cT.at(0).at(lengthOfSequence-1);
        for (int i=1; i < lengthOfSequence; i++){
            c.at(i).at(lengthOfSequence-1) = cT.at(i).at(lengthOfSequence-1) + c.at(i-1).at(lengthOfSequence-1);
        }
        // next fill the first row of c
        for (int j=lengthOfSequence-2; j>= 0; j--){
            c.at(0).at(j) = cT.at(0).at(j) + c.at(0).at(j+1);
            columnSumAbove.at(j) = cT.at(0).at(j);
        }
        // now fill the rest
        for (int i=1; i < lengthOfSequence; i++){
            for (int j=lengthOfSequence-2; j>= i; j--){
                c.at(i).at(j) = cT.at(i).at(j) + c.at(i).at(j+1) + columnSumAbove.at(j);
                columnSumAbove.at(j) += cT.at(i).at(j);
            }
        }

        return c;
    }

    // z-zeros
    vector<vector<int> > trimZeroOneZerosAllowed(
        const string& inputfile,
        const int& numberOfSequences,
        const int& lengthOfSequence,
        const int& numberOfAllowedZerosPerSequence,
        const int& thresholdGoodValues,
        const int& shiftToConvertChars)
    {
        int thresholdPlusShift = thresholdGoodValues + shiftToConvertChars;

        // alloc and init c and cC
        // cC = counter of columns
        // c(i,j)=m means, there are m lines where a block of "only ones and at most
        //          k zeros" starts at i and ends at j
        vector<vector<int> > c (lengthOfSequence, vector<int>(lengthOfSequence,0));
        vector<vector<int> > cC (lengthOfSequence, vector<int>(lengthOfSequence,0));

        // read the file row by row
        string zeile;

        // store the positions of the left ends of each 1-block
        vector<int> leftOne(lengthOfSequence,0);

        // store the positions of the right ends of each 1-block
        vector<int> rightOne(lengthOfSequence,0);

        // store the positions of all zeros in the current column
        vector<int> positionsOfZeros(lengthOfSequence,0);

        // if some zeros are allowed, there is a leftmost and a rightmost zero in
        // each identified block that consists of "only ones and at most k zeros"
        int leftBorderZero, rightBorderZero;
        // each block of such a type has a left and a right border
        int leftBorderOneBlock, rightBorderOneBlock;

        // open file
        ifstream in(inputfile, ios::in);

        bool stillInOneBlock;
        int startOfOneBlock;
        int numberOfZerosInCurrentRow;

        // loop over all lines of the file
        for (int z = 0; z < numberOfSequences; z++) {
            getline(in,zeile); // skip 3 lines
            getline(in,zeile); // skip 3 lines
            getline(in,zeile); // skip 3 lines
            getline(in,zeile);
            startOfOneBlock = 0;
            numberOfZerosInCurrentRow = 0;
            stillInOneBlock = (zeile.at(0)>=thresholdPlusShift);
            for (int i = 0; i < lengthOfSequence; i++) {
                // initialize leftOne and rightOne
                leftOne.at(i) = 0;
                rightOne.at(i) = 0;
                // store positions of Ones
                if (zeile.at(i)<thresholdPlusShift) {
                    positionsOfZeros.at(numberOfZerosInCurrentRow)=i;
                    numberOfZerosInCurrentRow++;
                }
                // fill leftOne
                if (zeile.at(i)<thresholdPlusShift) {
                    if (stillInOneBlock) {
                        stillInOneBlock = false;
                    }
                } else {
                    if (!stillInOneBlock) {
                        stillInOneBlock = true;
                        startOfOneBlock = i;
                    }
                    leftOne.at(i)=startOfOneBlock;
                }
            }
            // fill rightOne
            stillInOneBlock = (zeile.at(lengthOfSequence-1)>=thresholdPlusShift);
            if (stillInOneBlock) {
                startOfOneBlock = lengthOfSequence-1;
            }
            for (int i = lengthOfSequence-1; i >=0; i--) {
                if (zeile.at(i)<thresholdPlusShift) {
                    if (stillInOneBlock) {
                        stillInOneBlock = false;
                    }
                } else {
                    if (!stillInOneBlock) {
                        stillInOneBlock = true;
                        startOfOneBlock = i;
                    }
                    rightOne.at(i)=startOfOneBlock;
                }
            }

            if (numberOfZerosInCurrentRow <= numberOfAllowedZerosPerSequence) {
                for (int j=0; j < lengthOfSequence; j++)
                    cC.at(0).at(j)++;
            } else {
                int previousBlock = -1;
                for (int i = 0; i <= numberOfZerosInCurrentRow-numberOfAllowedZerosPerSequence; i++) {
                    leftBorderZero = positionsOfZeros.at(i);
                    rightBorderZero = positionsOfZeros.at(i+numberOfAllowedZerosPerSequence-1);
                    // leftBorderOneBlock is either
                    // 1) = 0, if leftBorderZero == 0
                    // 2) = leftBorderZero, if leftBorderZero-1 is *not* part of a
                    //                      1-block in zeile
                    // 3) = leftOne.at(leftBorderZero-1), if leftBorderZero-1 is part
                    //                                 of a 1-block in zeile
                    if ( leftBorderZero == 0 ) {
                        leftBorderOneBlock = 0;
                    } else { // leftBorderZero > 0
                        if (zeile.at(leftBorderZero-1)>=thresholdPlusShift) { // 1-block left of 0
                            leftBorderOneBlock = leftOne.at(leftBorderZero-1);
                        } else { // no 1-block
                            leftBorderOneBlock = leftBorderZero;
                        }
                    }
                    // same for rightBorderOneBlock
                    if (rightBorderZero == lengthOfSequence-1) {
                        rightBorderOneBlock = lengthOfSequence-1;
                    } else {
                        if (zeile.at(rightBorderZero+1)>=thresholdPlusShift) { // 1-block right of 0
                            rightBorderOneBlock = rightOne.at(rightBorderZero+1);
                        } else { // no 1-block
                            rightBorderOneBlock = rightBorderZero;
                        }
                    }
                    // add to cC
                    for (int j= previousBlock+1; j <=rightBorderOneBlock; j++) {
                        cC.at(leftBorderOneBlock).at(j)++;
                    }
                    previousBlock = rightBorderOneBlock;
                }
            }
        }

        // compute c from cC
        // first fill the first row of c
        for (int j=0; j< lengthOfSequence; j++){
            c.at(0).at(j) = cC.at(0).at(j);
        }
        // now fill the rest
        for (int i=1; i < lengthOfSequence; i++){
            for (int j=i; j<lengthOfSequence; j++){
                c.at(i).at(j) = cC.at(i).at(j) + c.at(i-1).at(j);
            }
        }

        return c;
    }

    // p-percent
    vector<vector<int> > trimZeroOnePercentZerosAllowed(
        const string& inputfile,
        const int& numberOfSequences,
        const int& lengthOfSequence,
        const double& percentOfAllowedZerosPerSequence,
        const int& thresholdGoodValues,
        const int& shiftToConvertChars)
    {
        int thresholdPlusShift = thresholdGoodValues + shiftToConvertChars;
        
        // c(i,j)=m means, there are m lines where a block of "only ones and at most
        //        p percent zeros" starts at i and ends at j
        vector<vector<int> > c (lengthOfSequence, vector<int>(lengthOfSequence,0));

        // read the file row by row
        string zeile;

        // pre compute allowed zeros per width for given percent
        vector<int> preCompAllowedZeros (lengthOfSequence+1);
        for (int i = 0; i <= lengthOfSequence; i++) {
            preCompAllowedZeros.at(i) = (int) (percentOfAllowedZerosPerSequence * i);
        }

        // row-block is feasible if numberOfCurrentZeros/(j-i+1) <= percentOfAllowedZerosPerSequence
        int numberOfCurrentZeros;
        int width;

        // open file
        ifstream in(inputfile, ios::in);

        for (int z = 0; z < numberOfSequences; z++) {
            getline(in,zeile); // skip 3 lines
            getline(in,zeile); // skip 3 lines
            getline(in,zeile); // skip 3 lines
            getline(in,zeile);
            for (int i = 0; i< lengthOfSequence; i++) {
                // init
                numberOfCurrentZeros = 0;
                width = 0;
                for (int j = i; j < lengthOfSequence; j++) {
                    width++;
                    if (zeile.at(j)<thresholdPlusShift) { numberOfCurrentZeros++; }
                    if (numberOfCurrentZeros <= preCompAllowedZeros.at(width) ) {
                        c.at(i).at(j)++;
                    }
                }
            }
        }

        return c;
    }

    // m-mean
    vector<vector<int> > trimIntegerMean(
        const string& inputfile,
        const int& numberOfSequences,
        const int& lengthOfSequence,
        const double& givenMean,
        const int& shiftToConvertChars)
    {
        // c(i,j)=x means, there are x lines where a block of starting at index i
        // and ending at index j with mean value at least "givenMean"
        vector<vector<int> > c (lengthOfSequence, vector<int>(lengthOfSequence,0));

        // read the file row by row
        string zeile;

        int currentSum;
        int currentWidth;
        double currentMean;
        double shiftedMean = shiftToConvertChars + givenMean;

        // open file
        ifstream in(inputfile, ios::in);

        for (int z = 0; z < numberOfSequences; z++) {
            getline(in,zeile); // skip 3 lines
            getline(in,zeile); // skip 3 lines
            getline(in,zeile); // skip 3 lines
            getline(in,zeile);
            for (int i = 0; i< lengthOfSequence; i++) {
                // init diagonal value
                currentSum = 0;
                currentWidth = 0;
                for (int j = i; j < lengthOfSequence; j++) {
                    currentSum += zeile.at(j);
                    currentWidth++;
                    currentMean = ( (double) currentSum ) / ( (double) currentWidth );
                    if ( currentMean >= shiftedMean ) { c.at(i).at(j)++; }
                }
            }
        }

        return c;
    }

}


/*******************************************************************************
 *
 * ComputeMatrices.h
 *
 * DESCRIPTION: Implementation of the algorithms for the problems:
 *              0-zeros:   trimZeroOne                    (line 43)
 *              z-zeros:   trimZeroOneZerosAllowed        (line 117)
 *              p-percent: trimZeroOnePercentZerosAllowed (line 269)
 *              m-mean:    trimIntegerMean                (line 402)
 *
 * RUNTIMES: If the input has r reads of length l:
 *           0-zeros:   worst-case: O( r * l )     expected: O( r * l )
 *           z-zeros:   worst-case: O( r * l )     expected: O( r * l )
 *           p-percent: worst-case: O( r * l^2 )   expected: O( r * l )
 *           m-mean:    worst-case: O( r * l^2 )   expected: O( r * l )
 *
 * AUTHORS: Ivo Hedtke (ivo.hedtke@uni-osnabrueck.de)
 *          Matthias Mueller-Hannemann (muellerh@informatik.uni-halle.de)
 *
 * CREATED: 13 Dec 2013
 *
 * LAST CHANGE: 25 Jul 2014
 *
 * Version 1.1: (25 Jul 2014) Speed-Up for p-percent and m-mean: If we have a
 *              nice instance, the runtime is O( r * l )
 *
 */

#include <fstream>
#include <string>
#include <vector>
#include <utility>

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
            in.ignore(numeric_limits<streamsize>::max(), '\n'); // skip 3 lines
            in.ignore(numeric_limits<streamsize>::max(), '\n'); // skip 3 lines
            in.ignore(numeric_limits<streamsize>::max(), '\n'); // skip 3 lines
            getline(in,zeile);
            startOfOneBlock = 0;
            zeile.append("!"); // dummy 0 at the end
            // ASSERT: ASCII("!")=33 and this is the smallest possible char in a quality score string
            // we need a dummy "bad" quality score at the end for our algorithm
            for (int i = 0; i <= lengthOfSequence; i++) {
                if (zeile[i] < thresholdPlusShift) {
                    if (stillInOneBlock) {
                        stillInOneBlock = false;
                        cT[startOfOneBlock][i-1]++;
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
        c[0][lengthOfSequence-1] = cT[0][lengthOfSequence-1];
        for (int i=1; i < lengthOfSequence; i++){
            c[i][lengthOfSequence-1] = cT[i][lengthOfSequence-1] + c[i-1][lengthOfSequence-1];
        }
        // next fill the first row of c
        for (int j=lengthOfSequence-2; j>= 0; j--){
            c[0][j] = cT[0][j] + c[0][j+1];
            columnSumAbove[j] = cT[0][j];
        }
        // now fill the rest
        for (int i=1; i < lengthOfSequence; i++){
            for (int j=lengthOfSequence-2; j>= i; j--){
                c[i][j] = cT[i][j] + c[i][j+1] + columnSumAbove[j];
                columnSumAbove[j] += cT[i][j];
            }
        }
        
        return c;
    }
    
    ////////////////////////////////////////////////////////////////////////////////
    
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
            in.ignore(numeric_limits<streamsize>::max(), '\n'); // skip 3 lines
            in.ignore(numeric_limits<streamsize>::max(), '\n'); // skip 3 lines
            in.ignore(numeric_limits<streamsize>::max(), '\n'); // skip 3 lines
            getline(in,zeile);
            startOfOneBlock = 0;
            numberOfZerosInCurrentRow = 0;
            stillInOneBlock = (zeile[0]>=thresholdPlusShift);
            for (int i = 0; i < lengthOfSequence; i++) {
                // initialize leftOne and rightOne
                leftOne[i] = 0;
                rightOne[i] = 0;
                // store positions of Ones
                if (zeile[i]<thresholdPlusShift) {
                    positionsOfZeros[numberOfZerosInCurrentRow]=i;
                    numberOfZerosInCurrentRow++;
                }
                // fill leftOne
                if (zeile[i]<thresholdPlusShift) {
                    if (stillInOneBlock) {
                        stillInOneBlock = false;
                    }
                } else {
                    if (!stillInOneBlock) {
                        stillInOneBlock = true;
                        startOfOneBlock = i;
                    }
                    leftOne[i]=startOfOneBlock;
                }
            }
            // fill rightOne
            stillInOneBlock = (zeile[lengthOfSequence-1]>=thresholdPlusShift);
            if (stillInOneBlock) {
                startOfOneBlock = lengthOfSequence-1;
            }
            for (int i = lengthOfSequence-1; i >=0; i--) {
                if (zeile[i]<thresholdPlusShift) {
                    if (stillInOneBlock) {
                        stillInOneBlock = false;
                    }
                } else {
                    if (!stillInOneBlock) {
                        stillInOneBlock = true;
                        startOfOneBlock = i;
                    }
                    rightOne[i]=startOfOneBlock;
                }
            }
            
            if (numberOfZerosInCurrentRow <= numberOfAllowedZerosPerSequence) {
                for (int j=0; j < lengthOfSequence; j++)
                    cC[0][j]++;
            } else {
                int previousBlock = -1;
                for (int i = 0; i <= numberOfZerosInCurrentRow-numberOfAllowedZerosPerSequence; i++) {
                    leftBorderZero = positionsOfZeros[i];
                    rightBorderZero = positionsOfZeros[i+numberOfAllowedZerosPerSequence-1];
                    // leftBorderOneBlock is either
                    // 1) = 0, if leftBorderZero == 0
                    // 2) = leftBorderZero, if leftBorderZero-1 is *not* part of a
                    //                      1-block in zeile
                    // 3) = leftOne[leftBorderZero-1], if leftBorderZero-1 is part
                    //                                 of a 1-block in zeile
                    if ( leftBorderZero == 0 ) {
                        leftBorderOneBlock = 0;
                    } else { // leftBorderZero > 0
                        if (zeile[leftBorderZero-1]>=thresholdPlusShift) { // 1-block left of 0
                            leftBorderOneBlock = leftOne[leftBorderZero-1];
                        } else { // no 1-block
                            leftBorderOneBlock = leftBorderZero;
                        }
                    }
                    // same for rightBorderOneBlock
                    if (rightBorderZero == lengthOfSequence-1) {
                        rightBorderOneBlock = lengthOfSequence-1;
                    } else {
                        if (zeile[rightBorderZero+1]>=thresholdPlusShift) { // 1-block right of 0
                            rightBorderOneBlock = rightOne[rightBorderZero+1];
                        } else { // no 1-block
                            rightBorderOneBlock = rightBorderZero;
                        }
                    }
                    // add to cC
                    for (int j= previousBlock+1; j <=rightBorderOneBlock; j++) {
                        cC[leftBorderOneBlock][j]++;
                    }
                    previousBlock = rightBorderOneBlock;
                }
            }
        }
        
        // compute c from cC
        // first fill the first row of c
        for (int j=0; j< lengthOfSequence; j++){
            c[0][j] = cC[0][j];
        }
        // now fill the rest
        for (int i=1; i < lengthOfSequence; i++){
            for (int j=i; j<lengthOfSequence; j++){
                c[i][j] = cC[i][j] + c[i-1][j];
            }
        }
        
        return c;
    }
    
    ////////////////////////////////////////////////////////////////////////////////
    
    //p-percent
    vector<vector<int> > trimZeroOnePercentZerosAllowed(const string& inputfile,
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
        vector<vector<int>> cT (lengthOfSequence, vector<int>(lengthOfSequence,0));
        
        // read the file row by row
        string zeile;
        
        // pre compute allowed zeros per width for given percent
        vector<int> preCompAllowedZeros (lengthOfSequence+1);
        for (int i = 0; i <= lengthOfSequence; i++) {
            preCompAllowedZeros[i] = (int) (percentOfAllowedZerosPerSequence * i);
        }
        
        // open file
        ifstream in(inputfile, ios::in);
        
        for (int z = 0; z < numberOfSequences; z++) {
            in.ignore(numeric_limits<streamsize>::max(), '\n'); // skip 3 lines
            in.ignore(numeric_limits<streamsize>::max(), '\n'); // skip 3 lines
            in.ignore(numeric_limits<streamsize>::max(), '\n'); // skip 3 lines
            getline(in,zeile);
            // pre processing to access the #zeros in O(1)
            // #zeros in g[L..R] equals partialSums[R+1] - partialSums[L]
            vector<int> partialSums(lengthOfSequence+1, 0);
            partialSums[0] = 0;
            for (int i = 0; i < lengthOfSequence; i++) {
                partialSums[i+1] = (zeile[i] < thresholdPlusShift) + partialSums[i];
            }
            
            // find block with values >= thresholdPlusShift), because all
            // subblocks fulfill the p-percent condition
            vector<pair<int,int>> oneBlocks;
            int startOfOneBlock = 0;
            bool stillInOneBlock = false;
            zeile.append("!"); // dummy 0 at the end
            // ASSERT: ASCII("!")=33 and this is the smallest possible char in a quality score string
            // we need a dummy "bad" quality score at the end for our algorithm
            for (int i = 0; i <= lengthOfSequence; i++) {
                if (zeile[i] < thresholdPlusShift) {
                    if (stillInOneBlock) {
                        stillInOneBlock = false;
                        oneBlocks.push_back( make_pair(startOfOneBlock,i-1) );
                        cT[startOfOneBlock][i-1]++;
                    }
                } else {
                    if (!stillInOneBlock) {
                        stillInOneBlock = true;
                        startOfOneBlock = i;
                    }
                }
            }
            
            // remember that we added "!" to zeile, but we will not read it later,
            // so there is no need to delete it
            
            // compute c(l,r) for all (l,r) not in the triangles of oneBlocks
            // HORIZONTAL
            int startrow = 0;
            for (auto p: oneBlocks) {
                for (int row = startrow; row < p.first; row++) {
                    for (int col = row+1; col < lengthOfSequence; col++) {
                        if ( (partialSums[col+1] - partialSums[row]) <= preCompAllowedZeros[col+1-row]) {
                            c[row][col]++;
                        }
                    }
                }
                startrow = p.second + 1;
                // VERTICAL: everything right of the triangle induced by p
                for (int row = p.first; row <= p.second; row++) {
                    for (int col = p.second+1; col < lengthOfSequence; col++) {
                        if ( (partialSums[col+1] - partialSums[row]) <= preCompAllowedZeros[col+1-row]) {
                            c[row][col]++;
                        }
                    }
                }
            }
            // everything after last triangle of 1s
            for (int row = startrow; row < lengthOfSequence; row++) {
                for (int col = row+1; col < lengthOfSequence; col++) {
                    if ( (partialSums[col+1] - partialSums[row]) <= preCompAllowedZeros[col+1-row]) {
                        c[row][col]++;
                    }
                }
            }
            
            
        }
        
        // compute c_aux from cT like in 0-zeros:
        vector<vector<int>> c_aux (lengthOfSequence, vector<int>(lengthOfSequence,0));
        vector<int> columnSumAbove (lengthOfSequence,0);
        // first fill the last column of c_aux
        c_aux[0][lengthOfSequence-1] = cT[0][lengthOfSequence-1];
        for (int i=1; i < lengthOfSequence; i++){
            c_aux[i][lengthOfSequence-1] = cT[i][lengthOfSequence-1] + c_aux[i-1][lengthOfSequence-1];
        }
        // next fill the first row of c_aux
        for (int j=lengthOfSequence-2; j>= 0; j--){
            c_aux[0][j] = cT[0][j] + c_aux[0][j+1];
            columnSumAbove[j] = cT[0][j];
        }
        // now fill the rest
        for (int i=1; i < lengthOfSequence; i++){
            for (int j=lengthOfSequence-2; j>= i; j--){
                c_aux[i][j] = cT[i][j] + c_aux[i][j+1] + columnSumAbove[j];
                columnSumAbove[j] += cT[i][j];
            }
        }
        
        // add c and c_aux
        
        for (int i = 0; i < lengthOfSequence; i++) {
            for (int j = i; j < lengthOfSequence; j++) {
                c[i][j] += c_aux[i][j];
            }
        }
        
        return c;
    }
    
    ////////////////////////////////////////////////////////////////////////////////
    
    // m-mean
    vector<vector<int>> trimIntegerMean(const string& inputfile,
                                        const int&    numberOfSequences,
                                        const int&    lengthOfSequence,
                                        const double& givenMean,
                                        const int&    shiftToConvertChars)
    {
        // c(i,j)=x means, there are x lines where a block of starting at index i
        // and ending at index j with mean value at least "givenMean"
        vector<vector<int>> c (lengthOfSequence, vector<int>(lengthOfSequence,0));
        vector<vector<int>> cT (lengthOfSequence, vector<int>(lengthOfSequence,0));
        
        double shiftedMean = shiftToConvertChars + givenMean;
        
        // read the file row by row
        string zeile;
        // open file
        ifstream in(inputfile, ios::in);
        
        for (int z = 0; z < numberOfSequences; z++) {
            in.ignore(numeric_limits<streamsize>::max(), '\n'); // skip 3 lines
            in.ignore(numeric_limits<streamsize>::max(), '\n'); // skip 3 lines
            in.ignore(numeric_limits<streamsize>::max(), '\n'); // skip 3 lines
            getline(in,zeile);
            
            // pre processing to access the mean in O(1)
            vector<int> partialSums(lengthOfSequence+1, 0);
            partialSums[0] = 0;
            for (int i = 0; i < lengthOfSequence; i++) {
                partialSums[i+1] = (zeile[i] - shiftedMean) + partialSums[i];
            }
            
            // find block with values >= mean, because all subblocks fulfill the
            // m-mean condition
            vector<pair<int,int>> oneBlocks;
            int startOfOneBlock = 0;
            bool stillInOneBlock = false;
            zeile.append("!"); // dummy 0 at the end
            // ASSERT: ASCII("!")=33 and this is the smallest possible char in a quality score string
            // we need a dummy "bad" quality score at the end for our algorithm
            for (int i = 0; i <= lengthOfSequence; i++) {
                if (zeile[i] < shiftedMean) {
                    if (stillInOneBlock) {
                        stillInOneBlock = false;
                        oneBlocks.push_back( make_pair(startOfOneBlock,i-1) );
                        cT[startOfOneBlock][i-1]++;
                    }
                } else {
                    if (!stillInOneBlock) {
                        stillInOneBlock = true;
                        startOfOneBlock = i;
                    }
                }
            }
            
            // remember that we added "!" to zeile, but we will not read it later,
            // so there is no need to delete it
            
            // compute c(l,r) for all (l,r) not in the triangles of oneBlocks
            // HORIZONTAL
            int startrow = 0;
            for (auto p: oneBlocks) {
                for (int row = startrow; row < p.first; row++) {
                    for (int col = row+1; col < lengthOfSequence; col++) {
                        if ( (partialSums[col+1] - partialSums[row]) >= 0) {
                            c[row][col]++;
                        }
                    }
                }
                startrow = p.second + 1;
                // VERTICAL: everything right of the triangle induced by p
                for (int row = p.first; row <= p.second; row++) {
                    for (int col = p.second+1; col < lengthOfSequence; col++) {
                        if ( (partialSums[col+1] - partialSums[row]) >= 0) {
                            c[row][col]++;
                        }
                    }
                }
            }
            // everything after last triangle of 1s
            for (int row = startrow; row < lengthOfSequence; row++) {
                for (int col = row+1; col < lengthOfSequence; col++) {
                    if ( (partialSums[col+1] - partialSums[row]) >= 0) {
                        c[row][col]++;
                    }
                }
            }
            
        }
        
        // compute c_aux from cT like in 0-zeros:
        vector<vector<int>> c_aux (lengthOfSequence, vector<int>(lengthOfSequence,0));
        vector<int> columnSumAbove (lengthOfSequence,0);
        // first fill the last column of c_aux
        c_aux[0][lengthOfSequence-1] = cT[0][lengthOfSequence-1];
        for (int i=1; i < lengthOfSequence; i++){
            c_aux[i][lengthOfSequence-1] = cT[i][lengthOfSequence-1] + c_aux[i-1][lengthOfSequence-1];
        }
        // next fill the first row of c_aux
        for (int j=lengthOfSequence-2; j>= 0; j--){
            c_aux[0][j] = cT[0][j] + c_aux[0][j+1];
            columnSumAbove[j] = cT[0][j];
        }
        // now fill the rest
        for (int i=1; i < lengthOfSequence; i++){
            for (int j=lengthOfSequence-2; j>= i; j--){
                c_aux[i][j] = cT[i][j] + c_aux[i][j+1] + columnSumAbove[j];
                columnSumAbove[j] += cT[i][j];
            }
        }
        
        // add c and c_aux
        
        for (int i = 0; i < lengthOfSequence; i++) {
            for (int j = i; j < lengthOfSequence; j++) {
                c[i][j] += c_aux[i][j];
            }
        }
        
        
        return c;
    }
    
}


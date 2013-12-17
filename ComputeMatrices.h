/*******************************************************************************
 *
 * ComputeMatrices.h
 *
 * AUTHORS: Ivo Hedtke (hedtke@informatik.uni-halle.de)
 *          Matthias Mueller-Hannemann (muellerh@informatik.uni-halle.de)
 *
 * CREATED: 13 Dec 2013
 *
 * LAST CHANGE: 16 Dec 2013
 *
 */

#include <fstream>

namespace ComputeMatrices {

    // free 2d array
    inline void free_matrix(int **matrix, int size_x) {
        for(int i = 0; i < size_x; i++) {
            free(matrix[i]);
        }
        free(matrix);
    }

    // calloc 1d array
    template<class T> inline T* callocVector(const int n) {
        return (T*) calloc(n,sizeof(T));
    }

    // malloc 1d array
    template<class T> inline T* mallocVector(const int n) {
        return (T*) malloc(n*sizeof(T));
    }

    // calloc 2d array
    template<class T> inline T** callocMatrix(const int n) {
        T** mat = mallocVector<T*>(n);
        for (int i = 0; i < n; i++) {
            mat[i] = callocVector<T>(n);
        }
        return mat;
    }

    // 0-zeros
    // =======
    // readFASTQ => each time skip 3 lines
    // (thresholdGoodValues >= 0) => lines of the grid are not '0' and '1' so use a threshold:
    //                               (quality < threshold)? '0' : '1'
    int** trimZeroOne(
            const std::string& inputfile,
            const int& numberOfSequences,
            const int& lengthOfSequence,
            const bool& readFASTQ,
            const int& thresholdGoodValues,
            const int& shiftToConvertChars)
    {
        // alloc and init c and cT
        // cT = counter of triangles
        // c(i,j)=m means, there are m lines where a 1-block starts at i & ends at j
        int** c = callocMatrix<int>(lengthOfSequence);
        int** cT = callocMatrix<int>(lengthOfSequence);

        // read the file row by row
        char* zeile = mallocVector<char>(lengthOfSequence+1);

        // open file
        std::ifstream in;
        in.open(inputfile.c_str(), std::ios::in);

        bool stillInOneBlock = false;
        int startOfOneBlock;
        for (int z = 0; z < numberOfSequences; z++) {
            if (readFASTQ) { // skip 3 lines
                in >> zeile;
                in >> zeile;
                in >> zeile;
            }
            in >> zeile;
            if (thresholdGoodValues >= 0) { // transform to 0/1
                for (int i = 0; i < lengthOfSequence; i++) {
                    zeile[i] = ( (zeile[i]-shiftToConvertChars)<thresholdGoodValues )? '0' : '1';
                }
            }
            startOfOneBlock = 0;
            zeile[lengthOfSequence] = '0'; // dummy 0 at the end
            for (int i = 0; i <= lengthOfSequence; i++) {
                if (zeile[i] == '0') {
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

        // close file
        in.close();

        // compute c from cT
        int* columnSumAbove = callocVector<int>(lengthOfSequence+1);
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

        free(zeile);
        free(columnSumAbove);
        free_matrix(cT, lengthOfSequence);

        return c;
    }

    // z-zeros
    int** trimZeroOneZerosAllowed(
            const std::string& inputfile,
            const int& numberOfSequences,
            const int& lengthOfSequence,
            const int& numberOfAllowedZerosPerSequence,
            const bool& readFASTQ,
            const int& thresholdGoodValues,
            const int& shiftToConvertChars)
    {

        // alloc and init c and cC
        // cC = counter of columns
        // c(i,j)=m means, there are m lines where a block of "only ones and at most
        //          k zeros" starts at i and ends at j
        int** c = callocMatrix<int>(lengthOfSequence);
        int** cC = callocMatrix<int>(lengthOfSequence);

        // read the file row by row
        char* zeile = mallocVector<char>(lengthOfSequence+1);

        // store the positions of the left ends of each 1-block
        int* leftOne = callocVector<int>(lengthOfSequence);

        // store the positions of the right ends of each 1-block
        int* rightOne = callocVector<int>(lengthOfSequence);

        // store the positions of all zeros in the current column
        int* positionsOfZeros = callocVector<int>(lengthOfSequence);

        // if some zeros are allowed, there is a leftmost and a rightmost zero in
        // each identified block that consists of "only ones and at most k zeros"
        int leftBorderZero, rightBorderZero;
        // each block of such a type has a left and a right border
        int leftBorderOneBlock, rightBorderOneBlock;

        // open file
        std::ifstream in;
        in.open(inputfile.c_str(), std::ios::in);

        bool stillInOneBlock;
        int startOfOneBlock;
        int numberOfZerosInCurrentRow;

        // loop over all lines of the file
        for (int z = 0; z < numberOfSequences; z++) {
            if (readFASTQ) { // skip 3 lines
                in >> zeile;
                in >> zeile;
                in >> zeile;
            }
            in >> zeile;
            if (thresholdGoodValues >= 0) { // transform to 0/1
                for (int i = 0; i < lengthOfSequence; i++) {
                    zeile[i] = ( (zeile[i]-shiftToConvertChars)<thresholdGoodValues )? '0' : '1';
                }
            }
            startOfOneBlock = 0;
            numberOfZerosInCurrentRow = 0;
            stillInOneBlock = (zeile[0]=='1');
            for (int i = 0; i < lengthOfSequence; i++) {
                // initialize leftOne and rightOne
                leftOne[i] = 0;
                rightOne[i] = 0;
                // store positions of Ones
                if (zeile[i]=='0') {
                    positionsOfZeros[numberOfZerosInCurrentRow]=i;
                    numberOfZerosInCurrentRow++;
                }
                // fill leftOne
                if (zeile[i] == '0') {
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
            stillInOneBlock = (zeile[lengthOfSequence-1]=='1');
            if (stillInOneBlock) {
                startOfOneBlock = lengthOfSequence-1;
            }
            for (int i = lengthOfSequence-1; i >=0; i--) {
                if (zeile[i] == '0') {
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
                        if (zeile[leftBorderZero-1]=='1') { // 1-block left of 0
                            leftBorderOneBlock = leftOne[leftBorderZero-1];
                        } else { // no 1-block
                            leftBorderOneBlock = leftBorderZero;
                        }
                    }
                    // same for rightBorderOneBlock
                    if (rightBorderZero == lengthOfSequence-1) {
                        rightBorderOneBlock = lengthOfSequence-1;
                    } else {
                        if (zeile[rightBorderZero+1]=='1') { // 1-block right of 0
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

        in.close();

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

        free(zeile);
        free(leftOne);
        free(rightOne);
        free(positionsOfZeros);
        free_matrix(cC, lengthOfSequence);

        return c;
    }

}
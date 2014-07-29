/*******************************************************************************
 *
 * ComputeMatricesParallel.h
 *
 * DESCRIPTION: Implementation of the algorithms for the problems:
 *              z-zeros:   trimZeroOneZerosAllowed
 *              p-percent: trimZeroOnePercentZerosAllowed
 *              m-mean:    trimIntegerMean
 *
 * RUNTIMES: If the input has r reads of length l:
 *           0-zeros:   worst-case: O( r * l )
 *           z-zeros:   worst-case: O( r * l )
 *           p-percent: worst-case: O( r * l^2 )
 *           m-mean:    worst-case: O( r * l^2 )
 *
 * AUTHORS: Ivo Hedtke (ivo.hedtke@uni-osnabrueck.de)
 *          Matthias Mueller-Hannemann (muellerh@informatik.uni-halle.de)
 *
 * CREATED: 25 Jul 2014
 *
 * LAST CHANGE: 29 Jul 2014
 *
 * Version 1.1: (29 Jul 2014) Parallel version
 *
 */

#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <thread>
#include <assert.h>

#include "ConcurrentQueue.h"

using namespace std;

namespace ComputeMatrices {
    
    void readFromFASTQFile(const string& inputfile,
                           const int& numberOfSequences,
                           ConcurrentQueue<string*>& q,
                           bool& ready){
        
        // this method reads a FASTQ-file line by line.
        // Every forth line (containing the quality information about a read)
        // is inserted into a thread-safe queue.
        // When parsing is completed, the booloean variable ready is set to true.
        
        // To keep extra space limited, the queue is limited to 1000 lines.
        // When the queue becomes full, the reading thread waits 2 miliseconds
        
        
        // read the file row by row
        string zeile;
        bool success = false;
        
        // open file
        ifstream in(inputfile.c_str(), ios::in);
        
        for (int z = 0; z < numberOfSequences; z++) {
            getline(in,zeile); // skip 3 lines
            getline(in,zeile); // skip 3 lines
            getline(in,zeile); // skip 3 lines
            getline(in,zeile);
            
            string* str = new string(zeile);
            
            success = false;
            
            while (!success){
                if (q.size() < 10000){
                    q.push(str);
                    success = true;
                }
                else {
                    //std::cout << "Queue is full\n";
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                }
            }
        }
        // parsing of the input file is completed
        ready = true;
        
    }

    /////////////////////////////////////////////////////////////////////////////
    // z-zeros
    
    void computeZeroOneZerosAllowedMatrix (ConcurrentQueue<string*>& q ,
                                                  vector<vector<int> >& c,
                                                  const int& lengthOfSequence,
                                                  const int& numberOfAllowedZerosPerSequence,
                                                  const int& thresholdPlusShift,
                                                  bool& ready)
    {

        // cC = counter of columns
        vector<vector<int> > cC (lengthOfSequence, vector<int>(lengthOfSequence,0));
        
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
      
        bool stillInOneBlock;
        int startOfOneBlock;
        int numberOfZerosInCurrentRow;
        
        while (!ready || !q.empty()){
            if (q.empty()){
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                //std::cout << "not ready, but queue is empty\n";
            }
            else {
                //if (!q.empty()){
                string* zeile = nullptr;
                
                if (!q.tryPop(zeile)){
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                    continue;
                }
                
                
                assert (zeile!= nullptr);
                
                //std::cout << "verarbeite zeile\n";
                
                startOfOneBlock = 0;
                numberOfZerosInCurrentRow = 0;
                stillInOneBlock = ((*zeile)[0]>=thresholdPlusShift);
                for (int i = 0; i < lengthOfSequence; i++) {
                    // initialize leftOne and rightOne
                    leftOne[i] = 0;
                    rightOne[i] = 0;
                    // store positions of Ones
                    if ((*zeile)[i]<thresholdPlusShift) {
                        positionsOfZeros[numberOfZerosInCurrentRow]=i;
                        numberOfZerosInCurrentRow++;
                    }
                    // fill leftOne
                    if ((*zeile)[i]<thresholdPlusShift) {
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
                stillInOneBlock = ((*zeile)[lengthOfSequence-1]>=thresholdPlusShift);
                if (stillInOneBlock) {
                    startOfOneBlock = lengthOfSequence-1;
                }
                for (int i = lengthOfSequence-1; i >=0; i--) {
                    if ((*zeile)[i]<thresholdPlusShift) {
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
                            if ((*zeile)[leftBorderZero-1]>=thresholdPlusShift) { // 1-block left of 0
                                leftBorderOneBlock = leftOne[leftBorderZero-1];
                            } else { // no 1-block
                                leftBorderOneBlock = leftBorderZero;
                            }
                        }
                        // same for rightBorderOneBlock
                        if (rightBorderZero == lengthOfSequence-1) {
                            rightBorderOneBlock = lengthOfSequence-1;
                        } else {
                            if ((*zeile)[rightBorderZero+1]>=thresholdPlusShift) { // 1-block right of 0
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
                
                delete zeile;
                
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
        
    }
    
    
    vector<vector<int> > trimZeroOneZerosAllowedPar(const string& inputfile,
                                                    const int& numberOfSequences,
                                                    const int& lengthOfSequence,
                                                    const int& numberOfAllowedZerosPerSequence,
                                                    const int& thresholdGoodValues,
                                                    const int& shiftToConvertChars,
                                                    const int& num_threads)
    {
        int thresholdPlusShift = thresholdGoodValues + shiftToConvertChars;
        
        // alloc and init c and cC
        // cC = counter of columns
        // c(i,j)=m means, there are m lines where a block of "only ones and at most
        //          k zeros" starts at i and ends at j
        vector<vector<int> > c (lengthOfSequence, vector<int>(lengthOfSequence,0));
 
        vector<vector <vector<int> > > cth (num_threads, vector< vector <int> >(lengthOfSequence, vector<int>(lengthOfSequence,0)));
        
        
        ConcurrentQueue<string*> q;
        bool ready = false;
        
        
        vector<thread> threads(num_threads);
        
        
        std::thread readerThread(std::bind(&readFromFASTQFile, inputfile, numberOfSequences, std::ref(q), std::ref(ready)));
        
        for (int i=0; i < num_threads; i++){
            threads[i] = thread(std::bind(&computeZeroOneZerosAllowedMatrix, std::ref(q), std::ref(cth[i]),lengthOfSequence,numberOfAllowedZerosPerSequence, thresholdPlusShift, std::ref(ready)));
        }
        
        // wait for all threads
        readerThread.join();
        std::for_each(threads.begin(), threads.end(),
                      std::mem_fn(&std::thread::join));
        
        // collect the results
        for (int th=0; th < num_threads; th++){
            for (int i = 0; i < lengthOfSequence; i++){
                for (int j=i; j <  lengthOfSequence; j++)
                    c[i][j] += cth[th][i][j];
            }
        }
        
        return c;
    
    }
    
    /////////////////////////////////////////////////////////////////////////////

    //p-percent
    
    void computeZeroOnePercentZerosAllowedMatrix (ConcurrentQueue<string*>& q ,
                                                  vector<vector<int> >& c,
                                                  const int& lengthOfSequence,
                                                  const double& percentOfAllowedZerosPerSequence,
                                                  const int& thresholdPlusShift,
                                                  bool& ready)
    {
 
        vector<vector<int>> cT (lengthOfSequence, vector<int>(lengthOfSequence,0));
        
        // pre compute allowed zeros per width for given percent
        vector<int> preCompAllowedZeros (lengthOfSequence+1);
        for (int i = 0; i <= lengthOfSequence; i++) {
            preCompAllowedZeros[i] = (int) (percentOfAllowedZerosPerSequence * i);
        }
        
        while (!ready || !q.empty()){
            if (q.empty()){
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                //std::cout << "not ready, but queue is empty\n";
            }
            else {
                //if (!q.empty()){
                string* zeile = nullptr;
                
                if (!q.tryPop(zeile)){
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                    continue;
                }
                
                
                assert (zeile!= nullptr);
                
                //std::cout << "verarbeite zeile\n";
 
                // pre processing to access the #zeros in O(1)
                // #zeros in g[L..R] equals partialSums[R+1] - partialSums[L]
                vector<int> partialSums(lengthOfSequence+1, 0);
                partialSums[0] = 0;
                for (int i = 0; i < lengthOfSequence; i++) {
                    partialSums[i+1] = ((*zeile)[i] < thresholdPlusShift) + partialSums[i];
                }
                
                // find block with values >= thresholdPlusShift), because all
                // subblocks fulfill the p-percent condition
                vector<pair<int,int>> oneBlocks;
                int startOfOneBlock = 0;
                bool stillInOneBlock = false;
                zeile->append("!"); // dummy 0 at the end
                // ASSERT: ASCII("!")=33 and this is the smallest possible char in a quality score string
                // we need a dummy "bad" quality score at the end for our algorithm
                for (int i = 0; i <= lengthOfSequence; i++) {
                    if ((*zeile)[i] < thresholdPlusShift) {
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

                delete zeile;

            }
        } // end while
 
        
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

        
    }
        
    
    vector<vector<int> > trimZeroOnePercentZerosAllowedPar(const string& inputfile,
                                                           const int& numberOfSequences,
                                                           const int& lengthOfSequence,
                                                           const double& percentOfAllowedZerosPerSequence,
                                                           const int& thresholdGoodValues,
                                                           const int& shiftToConvertChars,
                                                           const int num_threads )
    {
   
        int thresholdPlusShift = thresholdGoodValues + shiftToConvertChars;
        
        
        // c(i,j)=m means, there are m lines where a block of "only ones and at most
        //        p percent zeros" starts at i and ends at j
        vector<vector<int> > c (lengthOfSequence, vector<int>(lengthOfSequence,0));
        
        vector<vector <vector<int> > > cth (num_threads, vector< vector <int> >(lengthOfSequence, vector<int>(lengthOfSequence,0)));
        
        
        ConcurrentQueue<string*> q;
        bool ready = false;
        
        
        vector<thread> threads(num_threads);
        
        
        std::thread readerThread(std::bind(&readFromFASTQFile, inputfile, numberOfSequences, std::ref(q), std::ref(ready)));
        
        for (int i=0; i < num_threads; i++){
            threads[i] = thread(std::bind(&computeZeroOnePercentZerosAllowedMatrix, std::ref(q), std::ref(cth[i]),lengthOfSequence,percentOfAllowedZerosPerSequence, thresholdPlusShift, std::ref(ready)));
        }
        
        // wait for all threads
        readerThread.join();
        std::for_each(threads.begin(), threads.end(),
                      std::mem_fn(&std::thread::join));
        
        // collect the results
        for (int th=0; th < num_threads; th++){
            for (int i = 0; i < lengthOfSequence; i++){
                for (int j=i; j <  lengthOfSequence; j++)
                    c[i][j] += cth[th][i][j];
            }
        }
        
        return c;

        
    }
    
    /////////////////////////////////////////////////////////////////////////////
    
    void computeMeanMatrix (ConcurrentQueue<string*>& q ,
                            vector<vector<int> >& c,
                            const int& lengthOfSequence,
                            const double& givenMean,
                            const int& shiftToConvertChars,
                            bool& ready){
        
        int currentSum;
        double shiftedMean = shiftToConvertChars + givenMean;
        double currentCumulatedMean;
   
        vector<vector<int>> cT (lengthOfSequence, vector<int>(lengthOfSequence,0));
        
        // stop only if parsing is completed (ready == true) and
        // the queue has become empty (= every read has been processed)
        
        while (!ready || !q.empty()){
            if (q.empty()){
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                //std::cout << "not ready, but queue is empty\n";
            }
            else {
            //if (!q.empty()){
                string* zeile = nullptr;
                
                if (!q.tryPop(zeile)){
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                    continue;
                }
                
                
                assert (zeile!= nullptr);

                //std::cout << "verarbeite zeile\n";
                
                // pre processing to access the mean in O(1)
                vector<int> partialSums(lengthOfSequence+1, 0);
                partialSums[0] = 0;
                for (int i = 0; i < lengthOfSequence; i++) {
                    partialSums[i+1] = ((*zeile)[i] - shiftedMean) + partialSums[i];
                }
                
                // find block with values >= mean, because all subblocks fulfill the
                // m-mean condition
                vector<pair<int,int>> oneBlocks;
                int startOfOneBlock = 0;
                bool stillInOneBlock = false;
                zeile->append("!"); // dummy 0 at the end
                // ASSERT: ASCII("!")=33 and this is the smallest possible char in a quality score string
                // we need a dummy "bad" quality score at the end for our algorithm
                for (int i = 0; i <= lengthOfSequence; i++) {
                    if ((*zeile)[i] < shiftedMean) {
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
                
                delete zeile;
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
        
        
    }
   
    // m-mean (parallelized version)
    vector<vector<int> > trimIntegerMeanPar(
                                            const string& inputfile,
                                            const int& numberOfSequences,
                                            const int& lengthOfSequence,
                                            const double& givenMean,
                                            const int& shiftToConvertChars,
                                            const int num_threads   )
    {
        
        // c(i,j)=x means, there are x lines where a block of starting at index i
        // and ending at index j with mean value at least "givenMean"
        
        vector<vector<int> > c (lengthOfSequence, vector<int>(lengthOfSequence,0));
 
        
        vector<vector <vector<int> > > cth (num_threads, vector< vector <int> >(lengthOfSequence, vector<int>(lengthOfSequence,0)));
        
        
        ConcurrentQueue<string*> q;
        bool ready = false;
        
        
        vector<thread> threads(num_threads);
        
        
        std::thread readerThread(std::bind(&readFromFASTQFile, inputfile, numberOfSequences, std::ref(q),std::ref(ready)));
        
        for (int i=0; i < num_threads; i++){
            threads[i] = thread(std::bind(&computeMeanMatrix, std::ref(q), std::ref(cth[i]),lengthOfSequence,givenMean,shiftToConvertChars,std::ref(ready)));
        }
        
        // wait for all threads
        readerThread.join();
        std::for_each(threads.begin(), threads.end(),
                      std::mem_fn(&std::thread::join));
        
        // collect the results
        for (int th=0; th < num_threads; th++){
            for (int i = 0; i < lengthOfSequence; i++){
                for (int j=i; j <  lengthOfSequence; j++)
                    c[i][j] += cth[th][i][j];
            }
        }
        
        return c;
    }
    

    
    
}


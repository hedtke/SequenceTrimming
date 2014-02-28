/*******************************************************************************
 *
 * Results.h
 *
 * DESCRIPTION: Output routines. Given a matrix c with the number c(l,r) of
 *              reads that fulfill 0-zeros, z-zeros, p-percent or m-mean
 *              starting at column l and ending at column r.
 *              exportMatrix: exports the matrix c as CSV with the rows
 *                            "l; r; c(l,r)"
 *              printMaxArea: computes l',r':=argmax{ c(l,r)*(r-l+1) } and
 *                            returns c(l',r')*(r'-l'+1),
 *                                    r'-l'+1,
 *                                    c(l',r'),
 *                                    l',
 *                                    r'
 *                            on screen
 *
 * RUNTIMES: O(n^2) if the input matrix is of type (n x n).
 *
 * AUTHORS: Ivo Hedtke (ivo.hedtke@uni-osnabrueck.de)
 *          Matthias Mueller-Hannemann (muellerh@informatik.uni-halle.de)
 *
 * LAST CHANGE: 28 Feb 2014
 *
 */


using namespace std;

namespace Results {
    
    void exportMatrix(vector<vector<int> > c, string outfile) {
        ofstream out(outfile, ios::out);
        for ( int i = 0; i < c.size(); i++ ) {
            for ( int j = 0; j < c[i].size(); j++ ) {
                out << i << "; " << j << "; " c[i][j] << endl;
            }
        }
    }

    void printMaxArea(vector<vector<int> > c, int rows) {
        long long int maxvalue = 0;
        long long int value;
        int indexL = -1;
        int indexR = -1;
        for (int i = 0; i < c.size(); i++) {
            for (int j = i; j < c[i].size(); j++) {
                value = ((long long int) (j-i+1)) *((long long int) c[i][j]);
                if (value > maxvalue){
                    maxvalue = value;
                    indexL = i;
                    indexR = j;
                }
            }
        }
        int width = indexR - indexL + 1;
        cout << "area:  " << maxvalue << endl;
        cout << "width: " << width << " (" << (width*100.0)/((float) c.size()) << "%)" << endl;
        cout << "rows:  " << c[indexL][indexR] << " (" << (c[indexL][indexR]*100.0)/((float) rows) << "%)" << endl;
        cout << "left:  " << indexL << endl;
        cout << "right: " << indexR << endl;
    }

}
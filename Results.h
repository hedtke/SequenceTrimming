namespace Results {

    void exportResultMatrixAsCSV(int** c, int length, std::string outfile) {
        std::ofstream out;
        out.open(outfile.c_str(), std::ios::out);

        out << "left; right; rows;" << std::endl;
        for ( int i = 0; i < length; i++ ) {
            for ( int j = i; j < length; j++ ) {
                out << i << "; " << j << "; " << c[i][j] << ";" << std::endl;
            }
        }

        out.close();
    }

    void printMaxArea(int** c, int length, int rows) {
        long long int maxvalue = 0;
        long long int value;
        int indexL = -1;
        int indexR = -1;
        for (int i = 0; i < length; i++) {
            for (int j = i; j < length; j++) {
                value = (j-i+1) * c[i][j];
                if (value > maxvalue){
                    maxvalue = value;
                    indexL = i;
                    indexR = j;
                }
            }
        }
        int width = indexR - indexL + 1;
        std::cout << "area:  " << maxvalue << std::endl;
        std::cout << "width: " << width << " (" << (width*100.0)/((float) length) << "%)" << std::endl;
        std::cout << "rows:  " << c[indexL][indexR] << " (" << (c[indexL][indexR]*100.0)/((float) rows) << "%)" << std::endl;
        std::cout << "left:  " << indexL << std::endl;
        std::cout << "right: " << indexR << std::endl;
    }

}